import io, unittest
import vcf.model
import npsv.variant as variant
from npsv.feature_extraction import Variant

class VariantHelpersTestSuite(unittest.TestCase):
    """Variant helper functions"""

    def test_ci_for_precise_variants(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1 2827694 rs2376870 CGTGGATGCGGGGAC C . PASS SVTYPE=DEL;END=2827708;SVLEN=-14
"""
        )
        for record in vcf.Reader(vcf_file):
            self.assertTrue(record.is_sv)
            self.assertTrue(variant.is_precise(record))
            self.assertEqual(variant.get_ci(record, "CIPOS", 10), [0, 0])


    def test_ci_for_imprecise_variants(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
2 321682 . T <DEL> . PASS SVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20
"""
        )
        for record in vcf.Reader(vcf_file):
            self.assertTrue(record.is_sv)
            self.assertFalse(variant.is_precise(record))
            self.assertEqual(variant.get_ci(record, "CIPOS", 10), [-56, 20])
            self.assertEqual(variant.get_ci(record, "CIEND", 10), [-10, 10])


class VariantClassTestSuite(unittest.TestCase):
    def test_ci_for_precise_variants(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1 2827694 rs2376870 CGTGGATGCGGGGAC C . PASS SVTYPE=DEL;END=2827708;SVLEN=-14
"""
        )
        for record in vcf.Reader(vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record)
            self.assertTrue(variant.is_precise)
            self.assertEqual(variant.get_ci("CIPOS", 10), [0, 0])

    def test_ci_for_imprecise_variants(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
2 321682 . T <DEL> . PASS SVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20
"""
        )
        for record in vcf.Reader(vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record)
            self.assertFalse(variant.is_precise)
            self.assertEqual(variant.get_ci("CIPOS", 10), [-56, 20])
            self.assertEqual(variant.get_ci("CIEND", 10), [-10, 10])
