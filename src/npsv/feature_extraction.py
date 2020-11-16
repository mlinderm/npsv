import argparse, io, json, logging, os, re, subprocess, sys, tempfile
from collections import defaultdict
import vcf
import pysam
import pysam.bcftools as bcftools
import pybedtools.bedtool as bed
import numpy as np
import pandas as pd
from pathlib import Path
from shlex import quote
from .sample import Sample
from .variant import variant_descriptor, Variant, DeletionVariant
from npsv import npsva

ZSCORE_THRESHOLD = 1.5


def count_alleles_with_svviz2(variant, args, input_bam):
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
        "REF_READ",
        "ALT_READ",
        "REF_WEIGHTED_SPAN",
        "ALT_WEIGHTED_SPAN",
        "REF_CONC_SPAN",
        "ALT_CONC_SPAN",
        "INSERT_LOWER",
        "INSERT_UPPER",
        "CLIP_PRIMARY",
        "COVERAGE",
        "DHFC",
        "DHBFC",
        "DHFFC",
    ]
    FEATURE_COLS = ["#CHROM", "START", "END", "TYPE", "SAMPLE", "SVLEN", *FEATURES]
    MISSING_DATA = "."

    def __init__(self, variant: Variant, sample: Sample):
        self.__dict__['variant'] = variant
        self.__dict__['sample'] = sample
        self.__dict__['features'] = defaultdict(lambda: Features.MISSING_DATA)

    def __setattr__(self, name, value):
        if name not in Features.FEATURES:
            raise ValueError("Unknown feature: %s", name)
        self.features[name] = value

    def print_features(self, file, force_variant=None, ac=None):
        print_variant = self.variant if force_variant is None else force_variant
        print(
            print_variant.chrom,
            print_variant.pos,
            print_variant.end,
            print_variant.subtype,
            self.sample.name,
            print_variant.event_length,
            sep="\t",
            end="",
            file=file,
        )
        for feature in Features.FEATURES:
            file.write(f"\t{self.features[feature]}")
        if ac is not None:
            file.write(f"\t{ac}")
        file.write("\n")

    @staticmethod
    def header(out_file=sys.stdout, ac=None):
        """Generate header for SV features file

        Args:
            out_file (file object, optional): Output file object. Defaults to sys.stdout.
            ac (int, optional): Not None to include AC in header. Defaults to None.
        """
        if ac is not None:
            print(*Features.FEATURE_COLS, "AC", sep="\t", file=out_file)
        else:
            print(*Features.FEATURE_COLS, sep="\t", file=out_file)


def coverage_over_region(input_bam, region, min_mapq=40, min_baseq=15, min_anchor=11):
    """Compute coverage over 1-indexed, closed region"""
    depth_result = pysam.depth(  # pylint: disable=no-member
        "-Q", str(min_mapq),
        "-q", str(min_baseq),
        "-l", str(min_anchor),
        "-r", region,
        input_bam,
    )
    # start, end are 0-indexed half-open coordinates
    _, start, end = pysam.libcutils.parse_region(region=region)
    region_length = end - start
    if len(depth_result) > 0 and region_length > 0:
        depths = np.loadtxt(io.StringIO(depth_result), dtype=int, usecols=2)
        total_coverage = np.sum(depths)
        return (total_coverage / region_length, total_coverage, region_length)
    else:
        return (0., 0., region_length)

def count_realigned_reads(
    args,
    fragments,
    variant,
    ref_contig="ref",
    alt_contig="alt",
    count_straddle=True,
    **kwargs,
):

    # 1-indexed coordinates
    rl_breakpoint, rr_breakpoint = variant.ref_breakpoints(args.flank, contig=ref_contig)
    al_breakpoint, ar_breakpoint = variant.alt_breakpoints(args.flank, contig=alt_contig)

    counts, read_names = fragments.count_realigned_reads(
        [(rl_breakpoint, rr_breakpoint or "", al_breakpoint, ar_breakpoint or "")],
        count_straddle=count_straddle,
        **kwargs,
    )

    # If multiple breakpoints, average counts
    if args.mapq_reads:
        ref_reads = (counts["rl_mapq"] + counts["rr_mapq"]) / (1 if rr_breakpoint is None else 2)
        alt_reads = (counts["al_mapq"] + counts["ar_mapq"]) / (1 if ar_breakpoint is None else 2)
    else:
        ref_reads = (counts["rl"] + counts["rr"]) / (1 if rr_breakpoint is None else 2)
        alt_reads = (counts["al"] + counts["ar"]) / (1 if ar_breakpoint is None else 2)

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
        ref_contig=ref_contig,
        alt_contig=alt_contig,
        count_straddle=args.count_straddle,
    )
    features.REF_READ = ref_count
    features.ALT_READ = alt_count

    # Read Pair Evidence
    left_breakpoint = variant.left_flank_region_string(left_flank=1, right_flank=1)
    right_breakpoint = variant.right_flank_region_string(left_flank=1, right_flank=1)
    
    fragment_delta = -variant.event_length if variant.is_deletion else variant.event_length
    pair_results = fragments.count_pipeline_straddlers(
        left_breakpoint, right_breakpoint, pair_flank, fragment_delta, 1.5, args.min_anchor,
    )

    # TODO: Incorporate 'concordance' count features
    features.REF_WEIGHTED_SPAN = pair_results["ref_weighted_count"]
    features.ALT_WEIGHTED_SPAN = pair_results["alt_weighted_count"]
    features.REF_CONC_SPAN = pair_results["ref_conc_count"]
    features.ALT_CONC_SPAN = pair_results["alt_conc_count"]
    insert_count = pair_results["insert_count"]
    if insert_count > 0:
        features.INSERT_LOWER = pair_results["insert_lower"] / insert_count
        features.INSERT_UPPER = pair_results["insert_upper"] / insert_count
    else:
        features.INSERT_LOWER = 0
        features.INSERT_UPPER = 0

    # Clipped Read Evidence
    if variant.is_deletion:
        # For deletions we are interested in clipped reads within the event
        left_clip_results = fragments.count_pipeline_clipped_reads(left_breakpoint, args.min_clip)
        clip = left_clip_results["right"]
        clip_total = left_clip_results["total"]
        
        right_clip_results = fragments.count_pipeline_clipped_reads(right_breakpoint, args.min_clip)
        clip += right_clip_results["left"]
        clip_total += right_clip_results["total"]
    elif variant.is_insertion:
        # TODO: Handle complex variants
        # For insertions we are interested in clipped reads on either side of breakpoint
        clip_results = fragments.count_pipeline_clipped_reads(left_breakpoint, args.min_clip)
        clip = clip_results["left"] + clip_results["right"] + clip_results["both"]
        clip_total = clip_results["total"]
    elif variant.is_duplication:
        # For duplication we are interested in clipped reads outside the event
        left_clip_results = fragments.count_pipeline_clipped_reads(left_breakpoint, args.min_clip)
        clip = left_clip_results["left"]
        clip_total = left_clip_results["total"]
        
        right_clip_results = fragments.count_pipeline_clipped_reads(right_breakpoint, args.min_clip)
        clip += right_clip_results["right"]
        clip_total += right_clip_results["total"]

    features.CLIP_PRIMARY = clip / clip_total if clip_total > 0 else 0.


    # Coverage Evidence

    # Coverage currently only relevant for deletion events
    if not variant.is_deletion:
        return features

    coverage, _, _ = coverage_over_region(input_bam, variant.region_string(), min_mapq=args.min_mapq, min_baseq=args.min_baseq, min_anchor=args.min_anchor)
    _, left_flank_bases, left_flank_length = coverage_over_region(
        input_bam,
        variant.left_flank_region_string(args.rel_coverage_flank),
        min_mapq=args.min_mapq, min_baseq=args.min_baseq, min_anchor=args.min_anchor
    )
    _, right_flank_bases, right_flank_length = coverage_over_region(
        input_bam,
        variant.right_flank_region_string(args.rel_coverage_flank),
        min_mapq=args.min_mapq, min_baseq=args.min_baseq, min_anchor=args.min_anchor
    )

    # Normalized coverage features adapted from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6479422/

    # Chromosomal normalized coverage
    features.COVERAGE = coverage
    features.DHFC = coverage / sample.chrom_mean_coverage(variant.chrom)

    # GC normalized coverage
    gc_fraction, alignable_bases = variant.ref_gc_fraction
    if alignable_bases > 0:
        features.DHBFC = coverage / sample.gc_mean_coverage(gc_fraction)
    else:
        features.DHBFC = 1. if coverage > 0 else 0.

    # Flank normalized coverage
    total_flank_bases = left_flank_bases + right_flank_bases
    total_flank_length = left_flank_length + right_flank_length
    if total_flank_bases > 0 and total_flank_length > 0:
        features.DHFFC = coverage / (total_flank_bases / total_flank_length)
    else:
        features.DHFFC = 1. if coverage > 0 else 0.

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
        Features.header(out_file, ac)

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

