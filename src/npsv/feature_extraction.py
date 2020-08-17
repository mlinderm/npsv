import argparse, io, json, logging, os, re, subprocess, sys, tempfile
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
from .variant import variant_descriptor, Variant, DeletionVariant
from npsv import npsva

ZSCORE_THRESHOLD = 1.5


def count_alleles_with_npsva(
    args,
    variant,
    input_bam,
    sample: Sample,
    input_fasta=None,
    ref_contig="ref",
    alt_contig="alt",
    count_straddle=True,
    insert_hist=True,
    **kwargs,
):
    try:
        if input_fasta is None:
            # Generate FASTA with ref and alt alleles formatted for use with bwa (via npsva)
            fasta_path, ref_contig, alt_contig = variant.synth_fasta(args)
        else:
            fasta_path = input_fasta

        # Reference and alternate breakpoint spans in synthetic fasta (1-indexed)
        ref_length = variant.ref_length
        alt_length = variant.alt_length

        # 1-indexed coordinates
        rl_breakpoint = f"{ref_contig}:{args.flank}-{args.flank+1}"
        al_breakpoint = f"{alt_contig}:{args.flank}-{args.flank+1}"

        region = "{}:{}-{}".format(
            variant.record.CHROM,
            variant.record.POS - args.flank + 1,
            variant.end + args.flank,
        )

        count_alignment_args = {
            "region": region,
        }

        if ref_length > 1:
            count_alignment_args[
                "rr_region"
            ] = f"{ref_contig}:{args.flank + ref_length - 1}-{args.flank + ref_length}"

        if alt_length > 1:
            count_alignment_args[
                "ar_region"
            ] = f"{alt_contig}:{args.flank + alt_length - 1}-{args.flank + alt_length}"

        # Use or ignore actual density data (e.g. if simulated data)
        insert_size_density_dict = (
            sample.insert_size_density().as_dict() if insert_hist else {}
        )
        realigner = npsva.Realigner(
            fasta_path,
            sample.mean_insert_size,
            sample.std_insert_size,
            insert_size_density_dict,
        )

        counts, read_names = realigner.count_alignments(
            input_bam,
            rl_breakpoint,
            al_breakpoint,
            **count_alignment_args,
            count_straddle=count_straddle,
            **kwargs,
        )

        # If multiple alt breakpoints, average counts
        r_reads = (counts["rl_reads"] + counts["rr_reads"]) / (
            1 if ref_length == 1 else 2
        )
        a_reads = (counts["al_reads"] + counts["ar_reads"]) / (
            1 if alt_length == 1 else 2
        )

        return r_reads, a_reads, read_names
    finally:
        # Clean up the FASTA file if we created one
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

    def __init__(self, variant: Variant, sample: Sample):
        self.variant = variant
        self.sample = sample

    @property
    def read_counts(self):
        return (self.ref_reads, self.alt_reads)

    @read_counts.setter
    def read_counts(self, value):
        (self.ref_reads, self.alt_reads) = value

    def print_features(self, file, force_variant=None, ac=None):
        # TODO: Rationalize the column and attribute names
        print_variant = self.variant if force_variant is None else force_variant
        print(
            print_variant.chrom,
            print_variant.pos,
            print_variant.end,
            print_variant.subtype,
            self.sample.name,
            print_variant.event_length,
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



def extract_variant_features(
    args: argparse.Namespace,
    variant: Variant,
    input_bam: str,
    sample: Sample,
    max_reads: int = None,
    input_fasta: str = None,
    ref_contig: str = "ref",
    alt_contig: str = "alt",
    insert_hist: bool = True,
):
    with pysam.AlignmentFile(
        input_bam, mode="rb"
    ) as bam_reader:  # pylint: disable=no-member
        features = Features(variant, sample)

        # -------------------------------------
        # Split-read evidence
        # -------------------------------------
        ref_count, alt_count, *_ = count_alleles_with_npsva(
            args,
            variant,
            input_bam,
            sample,
            input_fasta=input_fasta,
            ref_contig=ref_contig,
            alt_contig=alt_contig,
            count_straddle=args.count_straddle,
            insert_hist=insert_hist,
        )
        # ref_count, alt_count, *_ = variant.count_alleles_with_svviz2(args, input_bam)
        features.read_counts = (ref_count, alt_count)

        # -------------------------------------
        # Paired-end evidence
        # -------------------------------------

        # TODO: Support insertions

        # Determine coordinates as 0-indexed [pos, end)
        # In correctly formatted VCF, POS is first base of event when zero-indexed, while
        # END is 1-indexed closed end or 0-indexed half-open end
        pos = variant.pos
        end = variant.end
        event_length = variant.ref_length

        ci_pos = variant.get_ci("CIPOS", default_ci=args.default_ci)
        ci_end = variant.get_ci("CIEND", default_ci=args.default_ci)

        # Determine paired-end evidence regions
        # TODO: Switch to sample.search_distance
        pair_flank = 1000

        left_id = right_id = bam_reader.gettid(variant.chrom)
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
                variant.chrom,
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
                variant.chrom,
                right_query_span,
                max_reads=max_reads,
            )
        else:
            query_span = (pos + ci_pos[0] - pair_flank, end + ci_pos[1] + pair_flank)
            gather_reads(
                fragments, bam_reader, variant.chrom, query_span, max_reads=max_reads
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
                print(fragment.query_name, ref_prob, alt_prob)
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
                print(fragment.query_name)
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

    
        # The following features are only relevant to deletions and duplications
        if not variant.is_deletion:
            return features

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
                f"{variant.chrom}:{start}-{stop}",
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
        coverage, _ = bases_in_region(variant.chrom, pos + 1, end)
        left_flank_coverage, _ = bases_in_region(
            variant.chrom, pos - args.rel_coverage_flank + 1, pos
        )
        right_flank_coverage, _ = bases_in_region(
            variant.chrom, end + 1, end + args.rel_coverage_flank
        )

        # Normalized coverage features adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6479422/

        # Chromosomal normalized coverage
        features.coverage = coverage
        features.dhfc = coverage / sample.chrom_mean_coverage(variant.chrom)

        # GC normalized coverage
        gc_fraction, alignable_bases = variant.ref_gc_fraction
        if alignable_bases > 0:
            features.dhbfc = coverage / sample.gc_mean_coverage(gc_fraction)
        else:
            # TODO: What does it imply for other coverage features that there are zero alignable bases?
            logging.warning("Zero alignable bases detected for variant %s", variant.id)
            features.dhbfc = 1.0

        # Flank normalized coverage
        mean_flank_coverage = (left_flank_coverage + right_flank_coverage) / 2
        if mean_flank_coverage > 0:
            features.dhffc = coverage / mean_flank_coverage
        else:
            features.dhffc = 1.0

        return features


def bases_in_region(input_bam, region, min_mapq=40, min_baseq=15, min_anchor=11):
    """Compute coverage over 1-indexed, closed region"""
    depth_result = pysam.depth(  # pylint: disable=no-member
        "-Q", str(min_mapq),
        "-q", str(min_baseq),
        "-l", str(min_anchor),
        "-r", region,
        input_bam,
    )
    depths = np.loadtxt(io.StringIO(depth_result), dtype=int, usecols=2)
    if depths.size > 0:
        _, start, end = pysam.libcutils.parse_region(region=region)
        mean_coverage = np.sum(depths) / (end - start)
        return mean_coverage
    else:
        return 0.

def count_realigned_reads(
    args,
    fragments,
    variant,
    ref_contig="ref",
    alt_contig="alt",
    count_straddle=True,
    **kwargs,
):

    # Reference and alternate breakpoint spans in synthetic fasta (1-indexed)
    ref_length = variant.ref_length
    alt_length = variant.alt_length

    # 1-indexed coordinates
    rl_breakpoint = f"{ref_contig}:{args.flank}-{args.flank+1}"
    rr_breakpoint = f"{ref_contig}:{args.flank + ref_length - 1}-{args.flank + ref_length}" if ref_length > 1 else ""
    al_breakpoint = f"{alt_contig}:{args.flank}-{args.flank+1}"
    ar_breakpoint = f"{alt_contig}:{args.flank + alt_length - 1}-{args.flank + alt_length}" if alt_length > 1 else ""

    counts, read_names = fragments.count_realigned_reads(
        [(rl_breakpoint, rr_breakpoint, al_breakpoint, ar_breakpoint)],
        count_straddle=count_straddle,
        **kwargs,
    )

    # If multiple breakpoints, average counts
    ref_reads = (counts["rl"] + counts["rr"]) / (1 if ref_length == 1 else 2)
    alt_reads = (counts["al"] + counts["ar"]) / (1 if alt_length == 1 else 2)

    return ref_reads, alt_reads, read_names


def extract_features(
    args: argparse.Namespace,
    variant: Variant,
    input_bam: str,
    sample: Sample,
    max_reads: int = None,
    input_fasta: str = None,
    ref_contig: str = "ref",
    alt_contig: str = "alt",
    insert_hist: bool = True,
):
    features = Features(variant, sample)

    try:
        if input_fasta is None:
            # Generate FASTA with ref and alt alleles formatted for use with bwa (via npsva)
            fasta_path, ref_contig, alt_contig = variant.synth_fasta(args)
        else:
            fasta_path = input_fasta

        fragments = npsva.RealignedFragments(
            fasta_path,
            sample.mean_insert_size,
            sample.std_insert_size,
            sample.insert_size_density().as_dict() if insert_hist else {},
            input_bam,
        )
    finally:
        # Clean up the FASTA file if we created one
        if fasta_path != input_fasta:
            os.remove(fasta_path)

    # Gather reads from the BAM file

    # Previously a fixed 1000 for paired-end and args.flank for realignment
    pair_flank = min(args.flank, sample.search_distance())  
    ci_pos = variant.get_ci("CIPOS", default_ci=args.default_ci)
    ci_end = variant.get_ci("CIEND", default_ci=args.default_ci)
    
    if variant.ref_length > 2 * pair_flank + ci_pos[1] + abs(ci_end[0]):
        # When event is large, gather reads only near the breakpoints
        fragments.gather_reads(variant.left_flank_region_string(left_flank=abs(ci_pos[0]) + pair_flank, right_flank=ci_pos[1] + pair_flank))       
        fragments.gather_reads(variant.right_flank_region_string(left_flank=abs(ci_end[0]) + pair_flank, right_flank=ci_end[1] + pair_flank))       
    else:
        # When event is small gather reads across the entire event
        fragments.gather_reads(variant.region_string(flank=pair_flank))        

    # Allele Count Evidence

    ref_count, alt_count, *_ = count_realigned_reads(
        args,
        fragments,
        variant,
        input_fasta=fasta_path,
        ref_contig=ref_contig,
        alt_contig=alt_contig,
        count_straddle=args.count_straddle,
    )
    features.read_counts = (ref_count, alt_count)

    # Read Pair Evidence
        
    pair_results = fragments.count_pipeline_straddlers(
        variant.region_string(), pair_flank, -variant.event_length, 1.5, args.min_anchor,
    )

    # TODO: Incorporate 'concordance' count features
    features.ref_span = pair_results["ref_weighted_count"]
    features.alt_span = pair_results["alt_weighted_count"]
    insert_count = pair_results["insert_count"]
    if insert_count > 0:
        features.insert_lower = pair_results["insert_lower"] / insert_count
        features.insert_upper = pair_results["insert_upper"] / insert_count
    else:
        features.insert_lower = 0
        features.insert_upper = 0

    # Coverage Evidence

    # Coverage currently only relevant for deletion events
    if not variant.is_deletion:
        return features

    coverage = bases_in_region(input_bam, variant.region_string(), min_mapq=args.min_mapq, min_baseq=args.min_baseq, min_anchor=args.min_anchor)
    left_flank_coverage = bases_in_region(
        input_bam,
        variant.left_flank_region_string(args.rel_coverage_flank),
        min_mapq=args.min_mapq, min_baseq=args.min_baseq, min_anchor=args.min_anchor
    )
    right_flank_coverage = bases_in_region(
        input_bam,
        variant.right_flank_region_string(args.rel_coverage_flank),
        min_mapq=args.min_mapq, min_baseq=args.min_baseq, min_anchor=args.min_anchor
    )

    # Normalized coverage features adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6479422/

    # Chromosomal normalized coverage
    features.coverage = coverage
    features.dhfc = coverage / sample.chrom_mean_coverage(variant.chrom)

    # GC normalized coverage
    gc_fraction, alignable_bases = variant.ref_gc_fraction
    if alignable_bases > 0:
        features.dhbfc = coverage / sample.gc_mean_coverage(gc_fraction)

    # Flank normalized coverage
    mean_flank_coverage = (left_flank_coverage + right_flank_coverage) / 2
    if mean_flank_coverage > 0:
        features.dhffc = coverage / mean_flank_coverage

    return features



def extract(
    args: argparse.Namespace,
    input_vcf: str,
    input_bam: str,
    out_file=sys.stdout,
    max_reads: int = None,
    ac: int = None,
    sample: Sample = None,
    force_variant: Variant = None,
    insert_hist: bool = True,
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
        force_variant (Variant, optional): Variant to determine feature variant columns. Defaults to None.
        insert_hist (bool, optional): Use empirical insert size historgram. Defaults to True.

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

    # Extract features for all SVs
    vcf_reader = vcf.Reader(filename=input_vcf)
    for record in vcf_reader:
        variant = Variant.from_pyvcf(record, args.reference)
        if variant is None:
            logging.warning("Variant type or VCF record not supported for %s. Skipping.", record.ID)

        features = extract_features(
            args,
            variant,
            input_bam,
            sample,
            max_reads=max_reads,
            insert_hist=insert_hist,
        )

        # Print features
        features.print_features(
            out_file, force_variant=force_variant, ac=ac,
        )

