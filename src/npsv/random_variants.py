import pandas as pd
import numpy as np
import pysam
import vcf
from operator import itemgetter

CHROM_REGEX_AUTO = r"^(chr)?([1-9][0-9]?)$"
CHROM_REGEX_AUTO_NOY = r"^(chr)?([1-9][0-9]?|[X])$"
CHROM_REGEX_ONLY_SEX = r"^(chr)?[XY]"


def sample_starts(contigs, size, n):
    linear_genome = np.zeros(contigs.shape[0] + 1, dtype=int)
    linear_genome[1:] = contigs["LENGTH"].cumsum()

    starts = np.random.random_integers(0, linear_genome[-1], size=n)
    starts.sort()
    bins = np.digitize(starts, linear_genome)

    start = starts - linear_genome[bins - 1]
    return pd.DataFrame(
        {
            "CHROM": contigs.CHROM.iloc[bins - 1],
            "START": start,
            "END": start + size,  # Add end coordinates (currently fixed size)
            "LINEAR": starts,
        },
        columns=["CHROM", "START", "END", "LINEAR"],
    )


def generate_n_variants(contigs, size, n, gaps):
    # Create variant regions, checking if they overlap known gaps
    variants = []
    while len(variants) < n:
        for variant in sample_starts(contigs, size, n=n - len(variants)).itertuples():
            # TODO: Check for variants beyond end of contig
            if variant[1] == 0:
                continue
            gap_iter = gaps.fetch(variant[1], variant[2], variant[3])
            if next(gap_iter, None):
                continue
            variants.append(variant[1:])  # Drop index field
    return variants


def random_variants(
    ref_path,
    genome_path,
    gap_path,
    output_file,
    size=300,
    n=1,
    use_X=False,
    only_sex=False,
    variant_path=None,
):
    # Filter contigs down as specified
    contigs = pd.read_table(
        genome_path, names=["CHROM", "LENGTH"], dtype={"CHROM": str}
    )
    filter_regex = CHROM_REGEX_AUTO
    if only_sex:
        filter_regex = CHROM_REGEX_ONLY_SEX
    elif use_X:
        filter_regex = CHROM_REGEX_AUTO_NOY
    contigs = contigs[contigs.CHROM.str.match(filter_regex)]

    gaps = pysam.TabixFile(gap_path)  # pylint: disable=no-member

    # If a VCF is specified, generate variants matching those in the VCF
    if variant_path is not None:
        variants = []
        for record in vcf.Reader(filename=variant_path):
            if not record.is_sv or record.var_subtype != "DEL":
                continue
            event_length = int(record.sv_end) - record.POS
            variants += generate_n_variants(contigs, event_length, 1, gaps)
    else:
        variants = generate_n_variants(contigs, size, n, gaps)

    # Sort variants in reference order (using linear coordinate)
    variants.sort(key=itemgetter(3))

    ref_reader = pysam.FastaFile(ref_path)  # pylint: disable=no-member

    print(
        """\
##fileformat=VCFv4.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">    
##ALT=<ID=DEL,Description="Deletion">""",
        "\n".join(
            map(
                lambda contig: "##contig=<ID={1},length={2}>".format(*contig),
                contigs.itertuples(),
            )
        ),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        sep="\n",
        file=output_file,
    )

    # Write VCF with variants
    for variant in variants:
        # Zero-indexed start is the same as the 1-indexed padding base (need to subtract one to fetch padding base)
        ref = ref_reader.fetch(variant[0], variant[1] - 1, variant[1])
        record = "{}\t{}\t.\t{}\t<DEL>\t.\t.\tSVTYPE=DEL;END={};SVLEN={};CIPOS=0,0;CIEND=0,0".format(
            variant[0], variant[1], ref, variant[2], -(variant[2] - variant[1])
        )
        print(record, file=output_file)

    ref_reader.close()
    gaps.close()

