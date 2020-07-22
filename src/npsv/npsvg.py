#!/usr/bin/env python3
import argparse, json, logging, os, re, sys, tempfile, vcf
from . import npsv_options


class HeaderAction(argparse.Action):
    def __init__(
        self,
        option_strings,
        dest=argparse.SUPPRESS,
        default=argparse.SUPPRESS,
        help="Print header line and exit",
    ):
        super(HeaderAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        from feature_extraction import header

        header(sys.stdout)
        parser.exit()


def main():
    parser = argparse.ArgumentParser(
        "npsvg", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Common options
    parser.add_argument(
        "-t",
        "--tempdir",
        help="Specify the temp directory",
        type=str,
        default=tempfile.gettempdir(),
    )
    parser.add_argument("--threads", help="Number of threads", type=int, default=1)
    logging_options = parser.add_mutually_exclusive_group()
    logging_options.add_argument(
        "-d",
        "--debug",
        help="Debug logging",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
        default=logging.WARNING,
    )
    logging_options.add_argument(
        "-v",
        "--verbose",
        help="Verbose logging",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
    )

    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")

    # Feature extraction
    parser_features = subparsers.add_parser("features", help="Extract features")
    parser_features.add_argument(
        "-i", "--input", help="Input VCF file.", type=str, dest="input", required=True
    )
    parser_features.add_argument(
        "-b", "--bam", help="Input BAM file.", type=str, dest="bam", required=True
    )
    parser_features.add_argument(
        "-o",
        "--output",
        action="store",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file",
    )
    parser_features.add_argument(
        "-r",
        "--reference-sequence",
        help="Reference fasta file.",
        type=str,
        dest="reference",
        required=True,
    )

    data_options = parser_features.add_argument_group()
    npsv_options.add_data_options(data_options)

    feature_options = parser_features.add_argument_group()
    npsv_options.add_feature_options(feature_options)
    feature_options.add_argument(
        "--no-header", dest="header", action="store_false", help="Suppress header line"
    )
    feature_options.add_argument("--header-only", action=HeaderAction)

    feature_options.add_argument("--ac", help="Set allele count", type=int)

    # Random variant generation
    parser_random = subparsers.add_parser("random", help="Extract features")
    parser_random.add_argument(
        "-i",
        "--input",
        help="Input VCF file to create equivalent random variants for",
        type=str,
        dest="input",
    )
    parser_random.add_argument(
        "-o",
        "--output",
        action="store",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file",
    )
    parser_random.add_argument(
        "-r",
        "--reference-sequence",
        help="Reference fasta file.",
        type=str,
        dest="reference",
        required=True,
    )
    parser_random.add_argument(
        "--genome", type=str, help="BedTools Genome file", required=True
    )
    npsv_options.add_random_options(parser_random)

    parser_random.add_argument(
        "--n", action="store", type=int, default=100, help="Number of variants"
    )

    # Genotyper
    parser_genotype = subparsers.add_parser("genotype", help="Genotype variants")
    parser_genotype.add_argument(
        "-i", "--input", help="Input VCF file.", type=str, dest="input", required=True
    )
    parser_genotype.add_argument(
        "--sim", type=str, help="Table of simulated data", required=True
    )
    parser_genotype.add_argument(
        "--real", type=str, help="Table of real data", required=True
    )
    parser_genotype.add_argument(
        "-o",
        "--output",
        action="store",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file",
    )
    parser_genotype.add_argument(
        "--samples",
        action="append",
        default=[],
        help="Force sample names when no calls in the file",
    )
    npsv_options.add_genotyping_options(parser_genotype)

    # Plotting
    parser_plot = subparsers.add_parser("plot", help="Plot features")
    parser_plot.add_argument(
        "--sim", type=str, help="Table of simulated data", required=True
    )
    parser_plot.add_argument(
        "--real", type=str, help="Table of real data", required=True
    )
    parser_plot.add_argument(
        "-v", "--vcf", type=str, help=".vcf file containing variants", required=True
    )
    parser_plot.add_argument(
        "-o", "--output", help="Output directory", type=str, required=True
    )

    # BAM preprocessing
    parser_preproc = subparsers.add_parser(
        "preprocess", help="Preprocess BAM to compute relevant stats"
    )
    parser_preproc.add_argument(
        "-b", "--bam", help="Input BAM file.", type=str, dest="bam", required=True
    )
    parser_preproc.add_argument(
        "-r",
        "--reference-sequence",
        help="Reference fasta file.",
        type=str,
        dest="reference",
        required=True,
    )
    parser_preproc.add_argument(
        "-o", "--output", type=str, help="Output file", required=True
    )
    parser_preproc.add_argument(
        "--goleft", default="goleft", type=str, help="Path to goleft executable",
    )
    parser_preproc_picard = parser_preproc.add_mutually_exclusive_group(required=True)
    parser_preproc_picard.add_argument(
        "--genome", type=str, help="BedTools Genome file"
    )
    parser_preproc_picard.add_argument(
        "--picard-gc",
        type=str,
        help="Path to Picard detail GC bias metrics",
        dest="picard_gc",
    )
    parser_preproc.add_argument(
        "--picard-insert",
        type=str,
        help="Path to Picard insert size metrics",
        dest="picard_insert",
    )
    parser_preproc.add_argument(
        "--picard-wgs", type=str, help="Path to Picard wgs metrics", dest="picard_wgs"
    )

    # Propose alternate representations
    parser_propose = subparsers.add_parser(
        "propose", help="Propose alternate representations"
    )
    parser_propose.add_argument(
        "-r",
        "--reference-sequence",
        help="Reference fasta file.",
        type=str,
        dest="reference",
        required=True,
    )
    parser_propose.add_argument(
        "-i", "--input", help="Input VCF file.", type=str, dest="input", required=True
    )
    parser_propose.add_argument(
        "-o",
        "--output",
        action="store",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file",
    )
    npsv_options.add_propose_options(parser_propose)

    # Select among proposed alternate representations
    parser_refine = subparsers.add_parser(
        "refine", help="Refine alternate representations to final VCF"
    )
    parser_refine.add_argument(
        "-i", "--input", help="Input VCF file.", type=str, dest="input", required=True
    )
    parser_refine.add_argument(
        "-o",
        "--output",
        action="store",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file",
    )

    # Generate FASTA with consensus sequence
    parser_consensus = subparsers.add_parser(
        "consensus", help="Generate FASTA with consensus sequence"
    )
    parser_consensus.add_argument(
        "-r",
        "--reference-sequence",
        help="Reference fasta file.",
        type=str,
        dest="reference",
        required=True,
    )
    parser_consensus.add_argument(
        "-i", "--input", help="Input VCF file.", type=str, dest="input", required=True
    )
    parser_consensus.add_argument(
        "-o",
        "--output",
        action="store",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file",
    )
    parser_consensus.add_argument(
        "--flank", help="Flank size for simulation region", type=int, default=3000
    )
    parser_consensus.add_argument(
        "--ac", help="Allele count for generating FASTA", type=int, default=1
    )
    parser_consensus.add_argument(
        "--ref_contig", help="Name for reference contig", type=str, default="ref"
    )
    parser_consensus.add_argument(
        "--alt_contig", help="Name for alternate contig", type=str, default="alt"
    )

    # GC-normalized simulation coverage
    parser_normcovg = subparsers.add_parser(
        "normcovg", help="Filter synthetic reads based on normalized GC coverage"
    )
    parser_normcovg.add_argument(
        "--fasta-path",
        dest="fasta_path",
        type=str,
        help="Path to FASTA file used to generate the synthetic reads",
        required=True
    )
    parser_normcovg.add_argument(
        "--stats-path",
        dest="stats_path",
        type=str,
        help="Path to stats JSON file generated by preprocessing command",
        required=True
    )
    parser_normcovg.add_argument(
        "-i", "--input", help="Input SAM file.", type=str, dest="input", required=True
    )
    npsv_options.add_simulation_options(parser_normcovg)

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=args.loglevel)

    if args.command == "features":
        from .feature_extraction import extract

        extract(args, args.input, args.bam, out_file=args.output, ac=args.ac, insert_hist=args.insert_hist)
    elif args.command == "genotype":
        from .genotyper import genotype_vcf

        genotype_vcf(
            args, args.input, args.sim, args.real, args.output, samples=args.samples
        )
    elif args.command == "random":
        from .random_variants import random_variants

        random_variants(
            args.reference,
            args.genome,
            args.gaps,
            args.output,
            size=args.size,
            n=args.n,
            use_X=args.use_X,
            only_sex=args.only_sex,
            variant_path=args.input,
        )
    elif args.command == "plot":
        from .plot import plot_features

        plot_features(args, args.sim, args.real, args.vcf, args.output)

    elif args.command == "preprocess":
        from .sample import compute_bam_stats

        # Compute stats and write to JSON file
        stats = compute_bam_stats(args, args.bam)
        with open(args.output, "w") as file:
            json.dump(stats, file)

    elif args.command == "propose":
        from .propose import propose_variants

        propose_variants(args, args.input, args.output)
    elif args.command == "refine":
        from .propose import refine_variants

        refine_variants(args, args.input, args.output)
    elif args.command == "consensus":
        from .variant import consensus_fasta

        consensus_fasta(args, args.input, args.output)

    elif args.command == "normcovg":
        from .simulation import filter_reads
        
        filter_reads(args, args.stats_path, args.fasta_path, args.input, "/dev/stdout")

if __name__ == "__main__":
    main()
