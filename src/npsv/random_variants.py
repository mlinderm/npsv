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


def generate_n_deletions(contigs, size, n, gaps, flank=0):
    # Create variant regions, checking if they overlap known gaps
    variants = []
    while len(variants) < n:
        for _, chrom, pos, end, linear in sample_starts(contigs, size, n=n - len(variants)).itertuples():
            assert chrom != 0, "Invalid chromosome"
            if pos < (2*flank + 1):
                continue
            if end >= (contigs[contigs["CHROM"] == chrom].iloc[0]["LENGTH"] - 2*flank):
                continue
            gap_iter = gaps.fetch(chrom, pos, end)
            if next(gap_iter, None):
                continue
            variants.append((chrom, pos, end, linear))
    return variants


def random_deletions(
    ref_path,
    genome_path,
    gap_path,
    output_file,
    size=300,
    n=1,
    use_X=False,
    only_sex=False,
    flank=0,
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

    
    # Generate and sort variants in reference order (using linear coordinate)
    try:
        variants = []
        for event_size in size:
            variants += generate_n_deletions(contigs, event_size, n, gaps, flank=flank)
    except TypeError:
        variants = generate_n_deletions(contigs, size, n, gaps, flank=flank)
    variants.sort(key=itemgetter(3))

    ref_reader = pysam.FastaFile(ref_path)  # pylint: disable=no-member

    # Write VCF with variants
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

    for variant in variants:
        # Zero-indexed start is the same as the 1-indexed padding base (need to subtract one to fetch padding base)
        ref = ref_reader.fetch(variant[0], variant[1] - 1, variant[1])
        record = "{}\t{}\t.\t{}\t<DEL>\t.\t.\tSVTYPE=DEL;END={};SVLEN={};CIPOS=0,0;CIEND=0,0".format(
            variant[0], variant[1], ref, variant[2], -(variant[2] - variant[1])
        )
        print(record, file=output_file)

    ref_reader.close()
    gaps.close()


def random_variants(
    variant_path,
    ref_path,
    genome_path,
    gap_path,
    output_file,
    n=1,
    use_X=False,
    only_sex=False,
    flank=0,
):
    sizes = []
    for record in vcf.Reader(filename=variant_path):
        if not record.is_sv or record.var_subtype != "DEL":
            continue
        event_length = int(record.sv_end) - record.POS
        sizes.append(event_length) 
    random_deletions(ref_path, genome_path, gap_path, output_file, size=sizes, n=n, use_X=use_X, only_sex=only_sex, flank=flank)
