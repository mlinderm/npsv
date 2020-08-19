import logging
import pandas as pd
import numpy as np
import pysam
import vcf
from operator import itemgetter
from npsv.variant import Variant

CHROM_REGEX_AUTO = r"^(chr)?([1-9][0-9]?)$"
CHROM_REGEX_AUTO_NOY = r"^(chr)?([1-9][0-9]?|[X])$"
CHROM_REGEX_SEX = r"^(chr)?[XY]"


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


def load_contigs(genome_path, use_X=False, only_sex=False):
    # Filter contigs down as specified
    contigs = pd.read_table(
        genome_path, names=["CHROM", "LENGTH"], dtype={"CHROM": str}
    )
    filter_regex = CHROM_REGEX_AUTO
    if only_sex:
        filter_regex = CHROM_REGEX_SEX
    elif use_X:
        filter_regex = CHROM_REGEX_AUTO_NOY
    contigs = contigs[contigs.CHROM.str.match(filter_regex)]
    return contigs    


def write_header(contigs, output_file):
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
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">""",
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


def generate_n_regions(contigs, size, n, gaps, flank=0):
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


def random_deletion(
    size,
    ref_reader, 
    contigs,
    gaps,
    output_file,
    n=1,
    flank=0,
):
    # Generate and sort variants in reference order (using linear coordinate)
    variants = generate_n_regions(contigs, size, n, gaps, flank=flank)
    variants.sort(key=itemgetter(3))

    for variant in variants:
        # Zero-indexed start of the deletion is the same as the 1-indexed padding base (need to subtract one to fetch padding base
        # when using 0-indexed interface to pysam)
        ref = ref_reader.fetch(variant[0], variant[1] - 1, variant[1])
        record = "{}\t{}\t.\t{}\t<DEL>\t.\t.\tSVTYPE=DEL;END={};SVLEN={};CIPOS=0,0;CIEND=0,0".format(
            variant[0], variant[1], ref, variant[2], -(variant[2] - variant[1])
        )
        print(record, file=output_file)


def random_insertion(
    insertion_seq,
    ref_reader, 
    contigs,
    gaps,
    output_file,
    n=1,
    flank=0,
):
    # Generate and sort variants in reference order (using linear coordinate)
    variants = generate_n_regions(contigs, 0, n, gaps, flank=flank)
    variants.sort(key=itemgetter(3))

    for variant in variants:
        # variant[1] is 1-indexed coordinate of padding base (need to subtract one to fetch padding base
        # when using 0-indexed interface to pysam)
        ref = ref_reader.fetch(variant[0], variant[1] - 1, variant[1])
        record = "{}\t{}\t.\t{}\t{}{}\t.\t.\tSVTYPE=INS;END={};SVLEN={};CIPOS=0,0;CIEND=0,0".format(
            variant[0], variant[1], ref, ref, insertion_seq, variant[1], len(insertion_seq)
        )
        print(record, file=output_file)


def random_variant(
    variant,
    ref_reader, 
    contigs,
    gaps,
    output_file,
    n=1,
    flank=0,
):
    # TODO: Handle complex variants more precisely, at present we only model the change in length (i.e. SVLEN)
    if variant.is_deletion:
        random_deletion(variant.event_length, ref_reader, contigs, gaps, output_file, n=n, flank=flank)
    elif variant.is_insertion:
        insertion_seq = variant._alt_seq(flank=0, ref_seq="")
        
        # insertion_seq should include the padding base
        assert len(insertion_seq) > variant.event_length
        random_insertion(insertion_seq[-variant.event_length:], ref_reader, contigs, gaps, output_file, n=n, flank=flank)


def random_variants(
    variant_or_vcf_path,
    ref_path,
    genome_path,
    gap_path,
    output_file,
    n=1,
    use_X=False,
    only_sex=False,
    flank=0,
):
    contigs = load_contigs(genome_path, use_X=use_X, only_sex=only_sex)
    gaps = pysam.TabixFile(gap_path)  # pylint: disable=no-member
    ref_reader = pysam.FastaFile(ref_path)  # pylint: disable=no-member

    write_header(contigs, output_file)

    if isinstance(variant_or_vcf_path, Variant):
        random_variant(variant_or_vcf_path, ref_reader, contigs, gaps, output_file, n=n, flank=flank)
    else:
        for record in vcf.Reader(filename=variant_or_vcf_path):
            variant = Variant.from_pyvcf(record)
            if variant is not None:
                random_variant(variant, ref_reader, contigs, gaps, output_file, n=n, flank=flank)
                
    ref_reader.close()
    gaps.close()