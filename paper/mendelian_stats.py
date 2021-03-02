#!/usr/bin/env python3
import argparse, sys
from enum import Enum
import vcf
import numpy as np
import pandas as pd
from npsv.variant import variant_descriptor

def make_argument_parser():
    """
    Construct argument parser
    """
    parser = argparse.ArgumentParser(
        __file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compute summary statistics from BCFTools Mendelian Errors",
    )
    parser.add_argument(
        "-t",
        help="Trio description (mother,father,proband)",
        type=str,
        dest="trio",
        required=True,
    )
    subparsers = parser.add_subparsers(dest="command", help="Sub-command help")
    parser_counts = subparsers.add_parser("counts", help="Count types of MERs")
    parser_gqs = subparsers.add_parser("gq", help="GQ of MERs")
    parser_list = subparsers.add_parser("list", help="MER sorted by minimum GQ")

    parser.add_argument("vcf", metavar="MER-VCF", type=str, help="BCFTools Mendelian error VCF file")
    return parser

# Types of Mendelian errors (proband, mother, father)
MER_GENOTYPES = {
    (1, 0, 0): "HetDenovo",
    (2, 0, 0): "HomAltDenovo",
    (0, 2, 0): "Undercall",
    (0, 2, 1): "Undercall",
    (0, 0, 2): "Undercall",
    (0, 1, 2): "Undercall",
    (0, 2, 2): "Undercall",
    (1, 2, 2): "Undercall",
    (2, 1, 0): "Overcall",
    (2, 2, 0): "Overcall",
    (2, 0, 1): "Overcall",
    (2, 0, 2): "Overcall",
}

def main(args):
    mother_sample, father_sample, proband_sample = args.trio.split(",")
    
    rows = []
    for record in vcf.Reader(filename=args.vcf):
        proband = record.genotype(proband_sample)
        mother = record.genotype(mother_sample)
        father = record.genotype(father_sample)
        
        genotypes = (proband.gt_type, mother.gt_type, father.gt_type)
        error_type = MER_GENOTYPES.get(genotypes, "Other")
        
        try:
            gqs = [call.data.GQ for call in (proband, mother, father) if "GQ" in call.data._fields]
            if len(gqs) == 0:
                gqs = [np.partition(call.data.PL, 1)[1] for call in (proband, mother, father) if "PL" in call.data._fields]
            if len(gqs) == 0:
                gqs = [call.data.SQ for call in (proband, mother, father) if "SQ" in call.data._fields]
            min_gq = np.min(gqs)
        except:
            min_gq = None

        rows.append((
            record.ID or variant_descriptor(record),
            error_type,
            min_gq,
        ))

    table = pd.DataFrame.from_records(
        rows, 
        columns=["ID","TYPE","MIN_GQ"],
    )
    table["TYPE"] = table["TYPE"].astype("category")           
    
    if args.command == "counts":
        type_counts = table.groupby("TYPE").size().to_frame("Count")
        type_counts["Fraction"] = type_counts["Count"] / type_counts["Count"].sum()
        type_counts.to_csv(sys.stdout)
    elif args.command == "gq":
        gq_bins = pd.cut(table.MIN_GQ, [0, 10, 20, 30, 40, 100], right=False)       
        gq_counts = table.groupby(gq_bins).size().to_frame("Count")
        gq_counts["Fraction"] = gq_counts["Count"] / gq_counts["Count"].sum()
        gq_counts.to_csv(sys.stdout, index_label="GQ")
    elif args.command == "list":
        list_table = table.sort_values(by=["MIN_GQ"], ascending=False).reset_index(drop=True)
        list_table["GIAB"] = False
        list_table.loc[(list_table["ID"] == "HG2_Ill_SpiralSDKrefine_6835") | (list_table["ID"] == "HG2_PB_SVrefine2PB10Xhap12_10613"),"GIAB"] = True
        list_table.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    parser = make_argument_parser()
    args = parser.parse_args()
    main(args)
