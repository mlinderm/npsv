import argparse


def add_random_options(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add arguments needed by random variant generator to parser
    
    Args:
        parser (argparse.ArgumentParser): Parser to extend to additional arguments
    
    Returns:
        argparse.ArgumentParser: Parser passed as argument
    """
    parser.add_argument(
        "--gaps",
        type=str,
        help="Tabix-indexed BED file of gaps in reference genome",
        required=True,
    )
    parser.add_argument(
        "--size", action="store", type=int, default=300, help="Size of variants"
    )
    parser.add_argument(
        "--use_X", action="store_true", default=False, help="Include X chromosome"
    )
    parser.add_argument(
        "--only_sex",
        action="store_true",
        default=False,
        help="Only generate variants on sex chromosomes",
    )
    return parser


def add_feature_options(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add arguments needed by feature extractor to parser
    
    Args:
        parser (argparse.ArgumentParser): Parser to extend to additional arguments
    
    Returns:
        argparse.ArgumentParser: Parser passed as argument
    """
    parser.add_argument(
        "--flank", help="Flank size for simulation region", type=int, default=3000
    )
    parser.add_argument(
        "--ci",
        dest="default_ci",
        default=10,
        type=int,
        help="Default CI for breakpoints",
    )
    parser.add_argument(
        "--min-anchor",
        dest="min_anchor",
        default=11,
        type=int,
        help="Min anchor span in bp",
    )
    parser.add_argument(
        "--min-mapQ",
        dest="min_mapq",
        default=40,
        type=int,
        help="Min read mapping quality to count reads towards depth",
    )
    parser.add_argument(
        "--min-baseQ",
        dest="min_baseq",
        default=15,
        type=int,
        help="Min base quality to count towards depth",
    )
    parser.add_argument(
        "--rel-covg-flank",
        dest="rel_coverage_flank",
        default=1000,
        type=int,
        help="Size of flanking regions for computing relative coverage to flanks",
    )
    return parser


def add_paragraph_options(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add arguments needed by paragraph to parser
    
    Args:
        parser (argparse.ArgumentParser): Parser to extend to additional arguments
    
    Returns:
        argparse.ArgumentParser: Parser passed as argument
    """
    parser.add_argument(
        "--vcf2paragraph",
        default="vcf2paragraph.py",
        type=str,
        help="Path to vcf2paragraph.py",
    )
    parser.add_argument("--grmpy", default="grmpy", type=str, help="Path to grmpy")
    parser.add_argument("--paragraph", default="paragraph", type=str, help="Path to paragraph")
    return parser


def add_data_options(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add arguments needed by simulation tools to parser
    
    Args:
        parser (argparse.ArgumentParser): Parser to extend to additional arguments
    
    Returns:
        argparse.ArgumentParser: Parser passed as argument
    """
    parser.add_argument(
        "-p", "--read-length", dest="read_length", type=int, help="Read length"
    )
    parser.add_argument(
        "--fragment-mean", dest="fragment_mean", type=float, help="Mean insert size"
    )
    parser.add_argument(
        "--fragment-sd",
        dest="fragment_sd",
        type=float,
        help="Standard deviation of fragment size",
    )
    parser.add_argument("--depth", type=float, help="Mean depth of coverage")
    parser.add_argument(
        "--stats-path",
        dest="stats_path",
        type=str,
        help="Path to stats JSON file generated by preprocessing command",
    )
    return parser


def add_genotyping_options(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add arguments needed by the genotyper to parser
    
    Args:
        parser (argparse.ArgumentParser): Parser to extend to additional arguments
    
    Returns:
        argparse.ArgumentParser: Parser passed as argument
    """
    parser.add_argument(
        "-c",
        "--classifier",
        help="Classifier(s) to use",
        type=str,
        default="svm"
    )
    parser.add_argument(
        "--local",
        help="Exclusively build per-variant classifier",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--downsample",
        type=int,
        help="Number of simulated replicates per variant and zygosity in 'global' mode",
        default=1,
    )
    parser.add_argument(
        "--filter-bed",
        dest="filter_bed",
        help="Only train on variants overlapping regions in this BED file",
        type=str,
    )
    parser.add_argument(
        "--dm2",
        help="Compute Mahalanobis distance for all genotypes",
        action="store_true",
        default=False,
    )
    return parser

def add_propose_options(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--simple-repeats-bed",
        dest="simple_repeats_bed",
        help="UCSC simple repeats BED file",
        type=str,
    )
    return parser