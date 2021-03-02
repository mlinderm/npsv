#!/usr/bin/env python3
import argparse, collections, sys
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
        "--sample",
        help="Sample name",
        type=str,
        default="HG002",
    )
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")
    parser_gq = subparsers.add_parser("gq", help="Concordance by GQ bucket")
    parser_vntr = subparsers.add_parser("vntr", help="Error enrichment by VNTR label")
    parser_len = subparsers.add_parser("len", help="Concordance by SVLEN")
    parser_list = subparsers.add_parser("list", help="List of all variants and metadata")
    parser_list = subparsers.add_parser("vntr-conc", help="Concordnace by VNTR label")
    parser.add_argument("vcf", metavar="TP-VCF", type=str,
                        help="Truvari 'tp-call.vcf' file")


    return parser


def main(args):
    # Construct pandas table from all VCF records
    rows = []
    regions = []
    for record in vcf.Reader(filename=args.vcf):
        call = record.genotype(args.sample)
        if "GQ" in call.data._fields:
            gq = call.data.GQ
        elif "PL" in call.data._fields and None not in call.data.PL:
            gq = np.partition(call.data.PL, 1)[1]
        elif "SQ" in call.data._fields:
            gq = call.data.SQ
        else:
            gq = None

        svlen = record.INFO["SVLEN"]
        if isinstance(svlen, collections.Sequence):
            svlen = svlen[0]

        rows.append((
            record.ID or variant_descriptor(record),
            record.INFO["SVTYPE"],
            abs(svlen),
            record.INFO["MatchGT"] == "TRUE",
            gq,
            record.INFO.get("TRall") == "TRUE",
            record.INFO.get("TRgt100") == "TRUE",
            record.INFO.get("TRgt10k") == "TRUE",
            record.INFO.get("segdup") == "TRUE",
            sum(map(lambda x: abs(x), record.INFO.get("CIPOS", [0,0]))),
            sum(map(lambda x: abs(x), record.INFO.get("CIEND", [0,0]))),
        ))

    table = pd.DataFrame.from_records(
        rows,
        columns=["ID", "TYPE", "SVLEN", "MatchGT", "GQ",
                 "TRall", "TRgt100", "TRgt10k", "SegDup", "CIPOS", "CIEND"],
    )

    if args.command == "gq":
        gq_bins = pd.cut(table.GQ, [0, 10, 20, 30, 40, 100], right=False)
        gq_table = table.groupby(gq_bins, observed=True).agg({"MatchGT": "value_counts"}).groupby(
            level=0).transform(lambda x: x / x.sum())

        idx = pd.IndexSlice
        gq_table.loc[idx[:, True], :].reset_index(level=1, drop=True).rename(
            columns={"MatchGT": "Concordance"}).to_csv(sys.stdout, index_label="GQ")
    elif args.command == "vntr":    
        enrichment = table.groupby("MatchGT").agg({
            "TRall": "value_counts",
            "TRgt100": "value_counts",
            "TRgt10k": "value_counts",
            "SegDup": "value_counts",
        }, normalize=True).dropna(axis=1, how="any").groupby(level=0).transform(lambda x: x/x.sum())
        
        idx = pd.IndexSlice
        enrichment.loc[idx[:, True], :].reset_index(level=1, drop=True).to_csv(
            sys.stdout, index_label=["Concordant"])
    elif args.command == "len":
        len_bins = pd.cut(
            table.SVLEN, [50, 100, 300, 1000, np.iinfo(np.int32).max], right=False)
        len_table = table.groupby(len_bins).agg({"MatchGT": "value_counts"}).groupby(
            level=0).transform(lambda x: x / x.sum())
        idx = pd.IndexSlice
        len_table.loc[idx[:, True], :].reset_index(level=1, drop=True).rename(
            columns={"MatchGT": "Concordance"}).to_csv(sys.stdout)
    elif args.command == "list":
        table.to_csv(sys.stdout, index=False)
    elif args.command == "vntr-conc":
        table["SVLEN_RANGE"] = pd.cut(table.SVLEN, [50, 100, 300, 1000, np.iinfo(np.int32).max], right=False)
        vntr_table = table.groupby(["SVLEN_RANGE","TRgt100"]).agg({"MatchGT": "value_counts"})
        
        vntr_table["Count"] = vntr_table.groupby(level=[0, 1]).transform(lambda x: x.sum())
        vntr_table["Concordance"] = vntr_table["MatchGT"] / vntr_table["Count"]

        idx = pd.IndexSlice
        vntr_table.loc[idx[:, :, True], :].reset_index(level=2, drop=True).drop(columns="MatchGT").to_csv(sys.stdout)


if __name__ == "__main__":
    parser = make_argument_parser()
    args = parser.parse_args()
    main(args)
