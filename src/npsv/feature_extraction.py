import argparse, json, logging, os, re, subprocess, sys, tempfile
import vcf
import pysam
import pysam.bcftools as bcftools
import pybedtools.bedtool as bed
import numpy as np
import pandas as pd
from pathlib import Path
from shlex import quote
from .sample import Sample
from .fragment import SpanningFragments, gather_reads
from .variant import get_ci, variant_descriptor, Variant, DeletionVariant
from npsv import npsva

ZSCORE_THRESHOLD = 1.5

def count_alleles_with_npsva(
    variant,
    args,
    input_bam,
    sample: Sample,
    input_fasta=None,
    ref_contig="ref",
    alt_contig="alt",
    **kwargs,
):
    try:
        if input_fasta is None:
            # Generate FASTA with ref and alt alleles formatted for use with bwa (via npsva)
            fasta_path, ref_contig, alt_contig = variant.synth_fasta(args)
        else:
            fasta_path = input_fasta

        # Reference and alternate breakpoint spans in synthetic fasta (1-indexed)
        length = variant.event_length
        alt_length = variant.alt_length

        rl_breakpoint = f"{ref_contig}:{args.flank}-{args.flank+1}"
        al_breakpoint = f"{alt_contig}:{args.flank}-{args.flank+1}"

        region = "{}:{}-{}".format(
            variant.record.CHROM,
            variant.record.POS - args.flank + 1,
            int(variant.record.sv_end) + args.flank,
        )

        count_alignment_args = {
            "rr_region": f"{ref_contig}:{args.flank + length}-{args.flank + length + 1}",
            "region": [region],
        }

        if alt_length > 1:
            count_alignment_args[
                "ar_region"
            ] = f"{alt_contig}:{args.flank+alt_length}-{args.flank+alt_length+1}"

        logging.debug(
            "Counting reads for %s and %s alleles in %s",
            ref_contig,
            alt_contig,
            fasta_path,
        )
        insert_size_density_dict = sample.insert_size_density().as_dict()
        realigner = npsva.Realigner(
            fasta_path,
            sample.mean_insert_size,
            sample.std_insert_size,
            insert_size_density_dict,
        )

        logging.debug(
            "Counting reads at ref. breakpoints (%s, %s) and alt. breakpoints (%s, %s)",
            rl_breakpoint,
            count_alignment_args["rr_region"],
            al_breakpoint,
            count_alignment_args.get("ar_region", None),
        )
        counts = realigner.count_alignments(
            input_bam,
            rl_breakpoint,
            al_breakpoint,
            **count_alignment_args,
            **kwargs,
        )

        # If multiple alt breakpoints, average counts
        r_reads = (counts["rl_reads"] + counts["rr_reads"]) / 2
        a_reads = (counts["al_reads"] + counts["ar_reads"]) / (
            1 if alt_length == 1 else 2
        )

        return r_reads, a_reads
    finally:
        # Clean up the file we created
        if fasta_path != input_fasta:
            os.remove(fasta_path)


def count_alleles_with_svviz2(variant: DeletionVariant, args, input_bam):
    with tempfile.TemporaryDirectory(dir=args.tempdir) as tempdir:
        vcf_path = variant.to_minimal_vcf(args, tempdir=tempdir)
        command = "{exec} --ref {ref} --outdir {outdir} --no-render --variants {vcf} {bam}".format(
            exec=quote("svviz2"),
            ref=quote(args.reference),
            outdir=tempdir,
            vcf=quote(vcf_path),
            bam=quote(input_bam),
        )
        subprocess.check_call(command, shell=True)

        id = (
            variant.record.ID
            or f"{variant.record.CHROM}_{variant.record.POS}_{variant.record.sv_end}"
        )
        svviz2_file_names = [
            f"{id}.{variant.record.var_subtype[:3].lower()}_{variant.record.CHROM}_{variant.record.POS-1}",
            f"{id}.{variant.record.var_subtype[:3].lower()}_{variant.record.CHROM}_{variant.record.POS}",
            f"{id}.SequenceDefinedVariant.{variant.record.CHROM}_{variant.record.POS-1}-{variant.record.sv_end-1}",
        ]
        for report_prefix in svviz2_file_names:
            report_path = os.path.join(tempdir, report_prefix + ".report.tsv")
            if os.path.exists(report_path):
                break
        else:
            assert False, f"Couldn't find report file for {id}"

        report = pd.read_csv(report_path, sep="\t")
        report.fillna({"allele": "all"}, inplace=True)

        alleles = report.groupby(["allele", "key"])
        ref = alleles.get_group(("ref", "count")).iat[0, 3]
        alt = alleles.get_group(("alt", "count")).iat[0, 3]

        return ref, alt


class Features(object):
    FEATURES = [
        "REF_SPAN",
        "ALT_SPAN",
        "REF_SPLIT",
        "ALT_SPLIT",
        "COVG",
        "INSERT_LOWER",
        "INSERT_UPPER",
        "DHFC",
        "DHBFC",
        "DHFFC",
    ]
    FEATURE_COLS = ["#CHROM", "START", "END", "TYPE", "SAMPLE", "SVLEN", *FEATURES]

    MISSING_DATA = "."

    def __init__(self, record: vcf.model._Record, sample: Sample):
        self.record = record
        self.sample = sample

    @property
    def read_counts(self):
        return (self.ref_reads, self.alt_reads)

    @read_counts.setter
    def read_counts(self, value):
        (self.ref_reads, self.alt_reads) = value

    def print_features(
        self, file, force_chrom=None, force_pos=None, force_end=None, ac=None
    ):
        # TODO: Rationalize the column and attribute names
        print(
            self.record.CHROM if force_chrom is None else force_chrom,
            self.record.POS if force_pos is None else force_pos,
            self.record.sv_end if force_end is None else force_end,
            self.record.var_subtype,
            self.sample.name,
            int(self.record.sv_end) - self.record.POS,
            getattr(self, "ref_span", Features.MISSING_DATA),
            getattr(self, "alt_span", Features.MISSING_DATA),
            getattr(self, "ref_reads", Features.MISSING_DATA),
            getattr(self, "alt_reads", Features.MISSING_DATA),
            getattr(self, "coverage", Features.MISSING_DATA),
            getattr(self, "insert_lower", Features.MISSING_DATA),
            getattr(self, "insert_upper", Features.MISSING_DATA),
            getattr(self, "dhfc", Features.MISSING_DATA),
            getattr(self, "dhbfc", Features.MISSING_DATA),
            getattr(self, "dhffc", Features.MISSING_DATA),
            sep="\t",
            end="",
            file=file,
        )
        if ac is not None:
            file.write(f"\t{ac}")
        file.write("\n")


def header(out_file=sys.stdout, ac=None):
    """Generate header for SV features file

    Args:
        out_file (file object, optional): Output file object. Defaults to sys.stdout.
        ac (int, optional): Not None to include AC in header. Defaults to None.
    """
    actual_features = Features.FEATURE_COLS[:]
    if ac is not None:
        actual_features.append("AC")
    print(*actual_features, sep="\t", file=out_file)


def prob_mapq(read):
    """Obtain mapping probability of read from MAPQ"""
    return 1 - 10 ** (-read.mapping_quality / 10.0)


def extract_id_from_paragraph(paragraph_result: dict) -> str:
    """Extract variant ID from sequence names in paragraph results

    Args:
        paragraph_result (dict): Paragraph output

    Returns:
        str: Variant ID
    """
    sequence_names = paragraph_result["sequencenames"][:]
    sequence_names.remove("REF")
    return os.path.commonprefix(sequence_names).rstrip(":")


def pad_vcf_alleles(args, input_vcf, output_vcf_file):
    vcf_reader = vcf.Reader(filename=input_vcf)
    records = [record for record in vcf_reader]
    assert len(records) == 1, "Can't pad more that one record"
    record = records[0]

    # Get the padding base from the reference
    fasta = pysam.FastaFile(args.reference)
    padding_base = fasta.fetch(record.CHROM, record.POS - 2, record.POS - 1)
    fasta.close()

    # Update the variant with a padding base
    record = vcf.model._Record(
        record.CHROM,
        record.POS - 1,
        record.ID,
        padding_base + record.REF,
        [vcf.model._Substitution(padding_base + str(alt)) for alt in record.ALT],
        record.QUAL,
        record.FILTER,
        record.INFO,
        record.FORMAT,
        record._sample_indexes,
        record.samples,
    )
    record_vcf_writer = vcf.Writer(output_vcf_file, vcf_reader)
    record_vcf_writer.write_record(record)


def convert_vcf_to_graph(args, input_vcf, output_graph):
    try:
        commandline = "{0} {1} -r {2} --graph-type alleles {3} {4}".format(
            quote(sys.executable),
            quote(args.vcf2paragraph),
            quote(args.reference),
            quote(input_vcf),
            quote(output_graph),
        )
        subprocess.check_call(commandline, shell=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        # The likely cause is overly strict parsing by Paragraph for non-simple insertions/deletions, add
        # in the desired padding base
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".vcf", dir=args.tempdir
        ) as padded_vcf_file:
            try:
                pad_vcf_alleles(args, input_vcf, padded_vcf_file)
                padded_vcf_file.close()

                commandline = "{0} {1} -r {2} --graph-type alleles {3} {4}".format(
                    quote(sys.executable),
                    quote(args.vcf2paragraph),
                    quote(args.reference),
                    quote(padded_vcf_file.name),
                    quote(output_graph),
                )

                subprocess.check_call(commandline, shell=True)
            finally:
                os.remove(padded_vcf_file.name)


def align_to_graph(
    args, input_bam, input_graphs, output_path, response_file: str = None
):
    if response_file is not None:
        resp_arg = f"--response-file {quote(response_file)}"
    else:
        resp_arg = ""

    commandline = "{exec} -r {ref} {graphs} -b {bam} -o {output} --log-level error --threads {threads} {resp} --variant-min-frac 0 --variant-min-reads 0".format(
        exec=quote(args.paragraph),
        ref=quote(args.reference),
        graphs=" ".join(map(lambda g: f"-g {quote(str(g))}", input_graphs)),
        bam=quote(input_bam),
        output=quote(output_path),
        threads=args.threads,
        resp=resp_arg,
    )
    subprocess.check_call(commandline, shell=True)


def run_paragraph_on_vcf(args, input_vcf: str, input_bam: str, sample: Sample):
    # Construct graph JSON file(s)
    if args.variant_json:
        converted_json_paths = [args.variant_json]
    else:
        converted_json_paths = []
        vcf_reader = vcf.Reader(filename=input_vcf)
        for record in vcf_reader:
            if not record.is_sv or record.var_subtype != "DEL":
                continue

            if record.ID is None:
                # Update ID to be more meaningful
                record.ID = variant_descriptor(record)

            try:
                # Split VCF into individual records
                with tempfile.NamedTemporaryFile(
                    mode="w", delete=False, suffix=".vcf", dir=args.tempdir
                ) as record_vcf_file:
                    record_vcf_writer = vcf.Writer(record_vcf_file.file, vcf_reader)
                    record_vcf_writer.write_record(record)

                # Generate graph JSON file for record
                record_json_path = Path(record_vcf_file.name).with_suffix(".json")
                convert_vcf_to_graph(args, record_vcf_file.name, str(record_json_path))
                converted_json_paths.append(record_json_path)
            finally:
                os.remove(record_vcf_file.name)

    # Run paragraph
    try:
        # Create paragraph response file with graph inputs to avoid issues with command line length
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".txt", dir=args.tempdir
        ) as response_file:
            for converted_json_path in converted_json_paths:
                print("-g", converted_json_path, sep=" ", file=response_file)

        paragraph_result_file = tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".json", dir=args.tempdir
        )
        paragraph_result_file.close()

        align_to_graph(
            args,
            input_bam,
            [],
            paragraph_result_file.name,
            response_file=response_file.name,
        )

        with open(paragraph_result_file.name, "r") as paragraph_result_file:
            results = json.load(paragraph_result_file)
            if not isinstance(results, list):
                results = [results]
            # Paragraph results aren't in a consistent order. Extract into dictionary keyed by variant ID.
            paragraph_dict = {}
            for result in results:
                variant_id = extract_id_from_paragraph(result)
                assert (
                    variant_id
                ), "Unable to extract valid variant ID from paragraph result"
                assert (
                    variant_id not in paragraph_dict
                ), "Duplicate variant ID in paragraph result"
                paragraph_dict[variant_id] = result
            return paragraph_dict
    finally:
        if not args.variant_json:
            map(os.remove, converted_json_paths)
        os.remove(response_file.name)
        os.remove(paragraph_result_file.name)


def extract_paragraph_read_counts(record, paragraph_output):
    id = record.ID
    assert id is not None, "Variant needs ID to robustly extract Paragraph results"
    ref_allele = id + ":0"
    alt_allele = id + ":1"

    read_counts = paragraph_output["read_counts_by_sequence"]
    reported_sequences = ref_sequence = alt_sequence = 0
    for allele, counts in read_counts.items():
        if ref_allele in allele:
            ref_sequence += counts["total:READS"]
            reported_sequences += 1
        if alt_allele in allele:
            alt_sequence += counts["total:READS"]
            reported_sequences += 1
    if reported_sequences == 0 and len(read_counts) > 0:
        logging.warning("No read counts detected for %s", variant_descriptor(record))

    return (ref_sequence, alt_sequence)


def extract(
    args: argparse.Namespace,
    input_vcf: str,
    input_bam: str,
    out_file=sys.stdout,
    max_reads: int = None,
    ac: int = None,
    sample: Sample = None,
    force_chrom: str = None,
    force_pos: int = None,
    force_end: int = None,
    input_fasta: str = None,
    ref_contig: str = "ref",
    alt_contig: str = "alt",
):
    """Extract and print deletion SV features for a VCF and BAM file

    Args:
        args (argparse.Namespace): Command line arguments
        input_vcf (str): Path to VCF file
        input_bam (str): Path to BAM file
        out_file (file_object, optional): File to write features. Defaults to sys.stdout.
        max_reads (int, optional): Max reads for feature extraction. Defaults to None.
        ac (int, optional): Allele count for current features. Defaults to None.
        sample (Sample, optional): Sample object. Defaults to None.
        force_chrom (str, optional): Overwrite actual CHROM if defined. Defaults to None.
        force_pos (int, optional): Overwrite actual POS if defined. Defaults to None.
        force_end (int, optional): Overwrite actual END if defined. Defaults to None.

    Raises:
        ValueError: Missing argument
    """
    # Print header line as requested
    if args.header:
        header(out_file, ac)

    if sample is not None:
        pass
    elif args.stats_path is not None:
        sample = Sample.from_npsv(args.stats_path, bam_path=input_bam)
    elif None not in (
        args.fragment_mean,
        args.fragment_sd,
        args.read_length,
        args.depth,
    ):
        sample = Sample.from_distribution(
            input_bam,
            args.fragment_mean,
            args.fragment_sd,
            args.read_length,
            mean_coverage=args.depth,
        )
    else:
        raise ValueError("Library distribution must be provided")

    # For efficiency we run paragraph on the entire VCF in one run (to amortize cost of loading the reference genome)
    # graph_alignments = run_paragraph_on_vcf(args, input_vcf, input_bam, sample)

    # Extract features for all SVs
    bam_reader = pysam.AlignmentFile(input_bam, mode="rb")  # pylint: disable=no-member
    vcf_reader = vcf.Reader(filename=input_vcf)
    for record in vcf_reader:
        if not record.is_sv:
            logging.warning("SVTYPE missing for variant %s. Skipping.", record.ID)
            continue

        kind = record.var_subtype
        if kind != "DEL":
            logging.warning("Unsupported SVTYPE for variant %s. Skipping", record.ID)
            continue

        # TODO: This code is in 3+ places, try to consolidate
        if record.ID is None:
            # Update ID to be more meaningful
            record.ID = variant_descriptor(record)

        features = Features(record, sample)

        # Determine coordinates as 0-indexed [pos, end)
        # In correctly formatted VCF, POS is first base of event when zero-indexed, while
        # END is 1-indexed closed end or 0-indexed half-open end
        pos = record.POS
        end = int(record.sv_end)
        event_length = end - pos

        ci_pos = get_ci(record, "CIPOS", default_ci=args.default_ci)
        ci_end = get_ci(record, "CIEND", default_ci=args.default_ci)

        # -------------------------------------
        # Split-read evidence
        # -------------------------------------
        variant = Variant.from_pyvcf(record)
        if variant is None:
            logging.warning("Unsupported variant type for %s. Skipping", record.ID)
            continue
        ref_count, alt_count, *_ = variant.count_alleles_with_npsva(
            args,
            input_bam,
            sample,
            input_fasta=input_fasta,
            ref_contig=ref_contig,
            alt_contig=alt_contig,
        )
        # ref_count, alt_count, *_ = variant.count_alleles_with_svviz2(args, input_bam)
        features.read_counts = (ref_count, alt_count)

        # assert record.ID in graph_alignments, "No paragraph data available"
        # ref_split, alt_split = extract_paragraph_read_counts(
        #     record, graph_alignments[record.ID]
        # )
        # features.read_counts = (ref_split, alt_split)

        # -------------------------------------
        # Paired-end evidence
        # -------------------------------------

        # Determine paired-end evidence regions
        pair_flank = 1000

        left_id = right_id = bam_reader.gettid(record.CHROM)
        left_paired_span = (pos + ci_pos[0] - pair_flank - 1, pos + ci_pos[1])
        right_paired_span = (end + ci_end[0], end + ci_pos[1] + pair_flank)

        fragments = SpanningFragments()
        if event_length > 2 * pair_flank + ci_pos[1] + abs(ci_end[0]):
            left_query_span = (
                pos + ci_pos[0] - pair_flank,
                pos + ci_pos[1] + pair_flank,
            )
            gather_reads(
                fragments,
                bam_reader,
                record.CHROM,
                left_query_span,
                max_reads=max_reads,
            )
            right_query_span = (
                end + ci_end[0] - pair_flank,
                end + ci_end[1] + pair_flank,
            )
            gather_reads(
                fragments,
                bam_reader,
                record.CHROM,
                right_query_span,
                max_reads=max_reads,
            )
        else:
            query_span = (pos + ci_pos[0] - pair_flank, end + ci_pos[1] + pair_flank)
            gather_reads(
                fragments, bam_reader, record.CHROM, query_span, max_reads=max_reads
            )

        # INSERT_UPPER and INSERT_LOWER adapted from:
        # https://www.sciencedirect.com/science/article/pii/S0092867418316337#sec4
        # ref_paired, alt_paired adapted from SVTyper
        ref_paired = alt_paired = insert_total = insert_upper = insert_lower = 0
        for fragment in fragments:

            # For deletion allow the concordant fragment to "jump" the entire event
            if fragment.is_pair_straddle(
                left_id, left_paired_span, right_id, right_paired_span, args.min_anchor
            ):
                alt_prob = fragment.prob_fragment_length(sample, event_length)
                ref_prob = fragment.prob_fragment_length(sample)
                if (alt_prob + ref_prob) != 0.0:
                    p_disc = alt_prob / (alt_prob + ref_prob)
                    alt_paired += p_disc
                    ref_paired += 1 - p_disc

                insert_total += 1
                zscore = fragment.zscore_fragment_length(sample)
                if zscore < -1.5:
                    insert_lower += 1
                elif zscore > 1.5:
                    insert_upper += 1

                continue

            # Right fragment of pair needs to be within the event
            left_ref_straddle = fragment.is_pair_straddle(
                left_id,
                (pos - pair_flank - 1, pos),
                left_id,
                (pos, min(end, pos + pair_flank)),
                args.min_anchor,
            )
            left_ref_straddle = (
                left_ref_straddle
                and pos + ci_pos[0] <= fragment.right.reference_start <= end + ci_end[1]
            )

            # Left fragment of pair needs to be within the event
            right_ref_straddle = fragment.is_pair_straddle(
                right_id,
                (max(pos, end - pair_flank - 1), end),
                right_id,
                (end, end + pair_flank),
                args.min_anchor,
            )
            right_ref_straddle = (
                right_ref_straddle
                and pos + ci_pos[0] <= fragment.left.reference_end <= end + ci_end[1]
            )

            if left_ref_straddle ^ right_ref_straddle:
                ref_paired += 0.5

                # Count all reference spanning reads towards insert_*
                insert_total += 1
                zscore = fragment.zscore_fragment_length(sample)
                if zscore < -1.5:
                    insert_lower += 1
                elif zscore > 1.5:
                    insert_upper += 1

        features.ref_span = ref_paired
        features.alt_span = alt_paired
        features.insert_lower = insert_lower / insert_total if insert_total != 0 else 0
        features.insert_upper = insert_upper / insert_total if insert_total != 0 else 0

        # -------------------------------------
        # Coverage evidence
        # -------------------------------------

        def bases_in_region(chrom, start, stop):
            """Compute coverage over 1-indexed, closed region [start,end]"""
            depth_result = pysam.depth(  # pylint: disable=no-member
                "-a",
                "-Q",
                str(args.min_mapq),
                "-q",
                str(args.min_baseq),
                "-l",
                str(args.min_anchor),
                "-r",
                f"{record.CHROM}:{start}-{stop}",
                bam_reader.filename,
            )

            depths = []
            for match in re.finditer(r"\S+\t\d+\t(\d+)\n", depth_result):
                if match:
                    depths.append(int(match.group(1)))
            site_length = stop - start + 1
            mean_coverage = np.sum(depths) / site_length if len(depths) > 0 else 0.0
            median_coverage = np.median(depths) if len(depths) > 0 else 0.0

            return mean_coverage, median_coverage

        # Use 1-indexed fully closed coordinates
        coverage, _ = bases_in_region(record.CHROM, pos + 1, end)
        left_flank_coverage, _ = bases_in_region(
            record.CHROM, pos - args.rel_coverage_flank + 1, pos
        )
        right_flank_coverage, _ = bases_in_region(
            record.CHROM, end + 1, end + args.rel_coverage_flank
        )

        # Normalized coverage features adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6479422/

        # Chromosomal normalized coverage
        features.coverage = coverage
        features.dhfc = coverage / sample.chrom_mean_coverage(record.CHROM)

        # GC normalized coverage

        # nucleotide_content columns: chrom, start, end, pct_at, pct_gc, num_A, num_C, num_G, num_T, num_N, num_oth, seq_len
        # pylint: disable=unexpected-keyword-arg
        nuc_content = bed.BedTool([(record.CHROM, pos, end)]).nucleotide_content(
            fi=args.reference
        )[0]
        alignable_bases = (
            int(nuc_content[11]) - int(nuc_content[10]) - int(nuc_content[9])
        )
        if alignable_bases > 0:
            gc_fraction = (int(nuc_content[6]) + int(nuc_content[7])) / alignable_bases
            features.dhbfc = coverage / sample.gc_mean_coverage(gc_fraction)
        else:
            # TODO: What does it imply for other coverage features that there are zero alignable bases?
            logging.warning("Zero alignable bases detected for variant %s", record.ID)
            features.dhbfc = 1.0

        # Flank normalized coverage
        logging.debug(
            "Coverage (left flank | event | right flank): %f | %f | %f",
            left_flank_coverage,
            coverage,
            right_flank_coverage,
        )
        mean_flank_coverage = (left_flank_coverage + right_flank_coverage) / 2
        if mean_flank_coverage > 0:
            features.dhffc = coverage / mean_flank_coverage
        else:
            features.dhffc = 1.0

        # -------------------------------------
        # Print features (including any derived features)
        # -------------------------------------
        features.print_features(
            out_file,
            force_chrom=force_chrom,
            force_pos=force_pos,
            force_end=force_end,
            ac=ac,
        )

    bam_reader.close()
