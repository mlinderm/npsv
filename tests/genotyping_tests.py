import argparse, io, os, sys, unittest
from npsv.genotyper import genotype_vcf

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class GenotypingCornerCasesTestSuite(unittest.TestCase):
    def test_empty_inputs(self):
        args = argparse.Namespace(local=False, filter_bed=None)
        out_file = io.StringIO()
        genotype_vcf(
            args,
            os.path.join(FILE_DIR, "empty.vcf"),
            os.path.join(FILE_DIR, "empty.tsv"),
            os.path.join(FILE_DIR, "empty.tsv"),
            output_file=out_file,
            samples=["HG003"],
        )
        header_line = out_file.getvalue().strip().split('\n')[-1]
        self.assertEqual(header_line, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003")
