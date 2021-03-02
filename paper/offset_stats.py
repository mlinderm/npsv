#!/usr/bin/env python3
import argparse, sys
import numpy as np
import pandas as pd
import vcf
import pybedtools.bedtool as bed
from npsv.variant import variant_descriptor

def make_argument_parser():
    """
    Construct argument parser
    """
    parser = argparse.ArgumentParser(
        __file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compute summary statistics from Truvari HG002 tp-call.vcf",
    )
    parser.add_argument(
        "-d",
        help="Offset table file",
        type=str,
        dest="offset_file",
        required=True,
    )
    
    parser.add_argument("vcf", metavar="TP-VCF", type=str, help="Truvari 'tp-call.vcf' file")

    return parser

def main(args):
    # Construct pandas table from all VCF records
    rows = []
    regions = []
    for record in vcf.Reader(filename=args.vcf):
        call = record.genotype("HG002")
        rows.append((
            record.ID or variant_descriptor(record),
            record.INFO["SVTYPE"],
            abs(record.INFO["SVLEN"][0]),
            record.INFO["MatchGT"] == "TRUE",
        ))

    table = pd.DataFrame.from_records(
        rows, 
        columns=["ID","TYPE","SIZE","MatchGT"],
    )
    
    # Potential offset (based on PBSV calls)
    if args.offset_file:
        offset_table = pd.read_table(args.offset_file)
        merge_table = table.merge(offset_table, on="ID")

        offset = merge_table[["StartDistance", "EndDistance"]].abs().max(axis=1)
        
        offset_intervals = pd.IntervalIndex.from_breaks([0, 1, 2, 5, 10, 20, 50, np.iinfo(np.int32).max], closed="left")
        offset_bins = pd.cut(offset, offset_intervals)
        offset_counts = merge_table.groupby(offset_bins).agg({ "MatchGT": "value_counts" }).groupby(level=0)
        
        totals = offset_counts.sum().rename(columns={ "MatchGT": "Total"})
        
        idx = pd.IndexSlice
        concordance = offset_counts.transform(lambda x: x / x.sum())
        concordance = concordance.loc[idx[:,True],:].reset_index(level=1, drop=True).rename(columns={"MatchGT": "Concordance"})
        
        totals.join(concordance).to_csv(sys.stdout, index_label="Offset")
   

if __name__ == "__main__":
    parser = make_argument_parser()
    args = parser.parse_args()
    main(args)
