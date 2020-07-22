import random, sys
import numpy as np
import pysam
from npsv.sample import Sample
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


def filter_reads(args, stats_path: str, fasta_path: str, in_sam: str, out_fastq: str):
    sample = Sample.from_npsv(stats_path)
    
    # Construct dictionary of GC normalized coverage
    library = sample.get_library(None)
    
    gc_dict = {}
    for gc in np.arange(0.0, 1.0, 0.01):
        gc_dict[round(gc * 100)] = library.gc_normalized_coverage(gc)

    max_normalized_gc = min(max(gc_dict.values()), args.max_gc_norm_covg)
    for gc, covg in gc_dict.items():
        gc_dict[gc] = covg / max_normalized_gc

    npsva.filter_reads(fasta_path, in_sam, out_fastq, gc_dict)

        