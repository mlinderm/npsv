import random, sys
import numpy as np
import pysam
import vcf
from npsv.sample import Sample
from npsv.variant import Variant
from npsv import npsva

def reverse_complement(sequence):
    result = ""
    for base in sequence[::-1]:
        if base == "A":
            result += "T"
        elif base == "C":
            result += "G"
        elif base == "G":
            result += "C"
        elif base == "T":
            result += "A"
    return result

def reverse(sequence, should_reverse):
    return sequence[::-1] if should_reverse else sequence


def filter_reads_gc(args, stats_path: str, fasta_path: str, in_sam: str, out_fastq: str):
    sample = Sample.from_npsv(stats_path)
    
    # Construct dictionary of GC normalized coverage
    library = sample.get_library(None)
    
    gc_dict = {}
    for gc in np.linspace(0.0, 1.0, 101):
        gc_dict[np.round(gc * 100).astype(int)] = library.gc_normalized_coverage(gc)
    max_normalized_gc = min(max(gc_dict.values()), args.max_gc_norm_covg)
    for gc, covg in gc_dict.items():
        gc_dict[gc] = covg / max_normalized_gc

    npsva.filter_reads_gc(fasta_path, in_sam, out_fastq, gc_dict)


def gnomad_coverage_profile(args, gnomad_coverage:str, input_vcf: str, output_file):
    record = next(vcf.Reader(filename=input_vcf))
    variant = Variant.from_pyvcf(record)
    variant.gnomad_coverage_profile(
        args,
        gnomad_coverage,
        output_file,
        ref_contig=args.ref_contig,
        alt_contig=args.alt_contig,
    )

# Actual max is 100, but for performance reasons we oversimulate a lower fraction
GNOMAD_MAX_COVG = 40. 
GNOMAD_MEAN_COVG = 30.6
GNOMAD_MULT_COVG = GNOMAD_MAX_COVG / GNOMAD_MEAN_COVG

def filter_reads_gnomad(args, covg_path: str, in_sam: str, out_fastq: str):
    npsva.filter_reads_gnomad(covg_path, in_sam, out_fastq, GNOMAD_MAX_COVG)
        