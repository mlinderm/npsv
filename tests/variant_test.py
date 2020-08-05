import argparse, io, os, sys, tempfile, unittest
from unittest.mock import patch
import vcf.model
from npsv.variant import Variant

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")

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
            variant = Variant.from_pyvcf(record, None)
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
            variant = Variant.from_pyvcf(record, None)
            self.assertFalse(variant.is_precise)
            self.assertEqual(variant.get_ci("CIPOS", 10), [-56, 20])
            self.assertEqual(variant.get_ci("CIEND", 10), [-10, 10])


class SimpleVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1	899922	HG3_Ill_SVrefine2DISCOVARDovetail_2	GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGT	G	20	PASS	END=899992;SVTYPE=DEL;SVLEN=-70
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            self.assertEqual(ref_contig, "1_899922_899993")
            self.assertEqual(alt_contig, "1_899922_899993_alt")
            mock_ref.assert_called_once_with(region="1:899922-899993")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">1_899922_899993")
            self.assertEqual(
                lines[1],
                "GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
            )
            self.assertEqual(lines[2], ">1_899922_899993_alt")
            self.assertEqual(lines[3], "GG")

    def test_consensus_fasta_contig_names(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(
                self.args, ref_contig="ref", alt_contig="alt", line_width=sys.maxsize
            )
            self.assertEqual(ref_contig, "ref")
            self.assertEqual(alt_contig, "alt")
            mock_ref.assert_called_once_with(region="1:899922-899993")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">ref")
            self.assertEqual(lines[2], ">alt")

    def test_consensus_fasta_hom_ref(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(
                self.args, ac=0, ref_contig="ref", alt_contig="alt", line_width=sys.maxsize
            )
            self.assertEqual(ref_contig, "ref")
            self.assertEqual(alt_contig, "alt")
            mock_ref.assert_called_once_with(region="1:899922-899993")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 3)
            self.assertEqual(lines[0], ">ref")
            self.assertEqual(
                lines[1],
                "GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
            )
            self.assertEqual(lines[2], ">alt")

    def test_consensus_fasta_hom_alt(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(
                self.args, ac=2, ref_contig="ref", alt_contig="alt", line_width=sys.maxsize
            )
            self.assertEqual(ref_contig, "ref")
            self.assertEqual(alt_contig, "alt")
            mock_ref.assert_called_once_with(region="1:899922-899993")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 3)
            self.assertEqual(lines[0], ">ref")
            self.assertEqual(lines[1], ">alt")
            self.assertEqual(lines[2], "GG")

    def test_gnomad_coverage_profile(self):
        record = next(vcf.Reader(self.vcf_file))
        self.assertTrue(record.is_sv)
        variant = Variant.from_pyvcf(record, None)

        fasta_path, ref_contig, alt_contig = variant.gnomad_coverage_profile(
            self.args, 
            os.path.join(FILE_DIR, "1_896922_903086.gnomad.genomes.coverage.summary.tsv.gz"),
            line_width=sys.maxsize,
        )
        self.assertEqual(ref_contig, "1_899922_899993")
        self.assertEqual(alt_contig, "1_899922_899993_alt")
        with open(fasta_path, "r") as fasta:
            lines = [line.strip() for line in fasta]
        self.assertEqual(len(lines), 4)
        self.assertEqual(lines[0], ">1_899922_899993")
        self.assertEqual(lines[2], ">1_899922_899993_alt")
        self.assertEqual(lines[3], "=1")


class SymbolicVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1	899922	HG3_Ill_SVrefine2DISCOVARDovetail_2	G	<DEL>	20	PASS	END=899992;SVTYPE=DEL;SVLEN=-70
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            self.assertEqual(ref_contig, "1_899922_899993")
            self.assertEqual(alt_contig, "1_899922_899993_alt")
            mock_ref.assert_called_once_with(region="1:899922-899993")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">1_899922_899993")
            self.assertEqual(
                lines[1],
                "GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGTG",
            )
            self.assertEqual(lines[2], ">1_899922_899993_alt")
            self.assertEqual(lines[3], "GG")


class ComplexVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
4	20473846	HG3_PB_assemblyticsfalcon_19066	TATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATA    ATATATATAGATATATCTATATATATCTATATAGATATATCTATATCTATATAGATATATCTATATATATATAGATATATCTATATCTATATAGATATATATCTATATATATATCTATATAGATATATCTATATAGATATAGATATATATCTATATATAGATATAGATATATCTATATAGATATATATCTATAGATATCTATATATATAGATATATAGATATCTATATCTATAT	10   LongReadHomRef	END=20474269;SVTYPE=DEL;SVLEN=-190
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        self.assertTrue(record.is_sv)
        variant = Variant.from_pyvcf(record, None)
        self.assertTrue(variant.is_deletion)
        self.assertEqual(
            variant.alt_length,
            len(
                "ATATATATAGATATATCTATATATATCTATATAGATATATCTATATCTATATAGATATATCTATATATATATAGATATATCTATATCTATATAGATATATATCTATATATATATCTATATAGATATATCTATATAGATATAGATATATATCTATATATAGATATAGATATATCTATATAGATATATATCTATAGATATCTATATATATAGATATATAGATATCTATATCTATAT"
            ),
        )

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="TATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATAT",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            mock_ref.assert_called_once_with(region="4:20473846-20474270")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">4_20473846_20474270")
            self.assertEqual(lines[0], f">{ref_contig}")
            self.assertEqual(
                lines[1],
                "TATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATAT",
            )
            self.assertEqual(lines[2], ">4_20473846_20474270_alt")
            self.assertEqual(lines[2], f">{alt_contig}")
            self.assertEqual(
                lines[3],
                "ATATATATAGATATATCTATATATATCTATATAGATATATCTATATCTATATAGATATATCTATATATATATAGATATATCTATATCTATATAGATATATATCTATATATATATCTATATAGATATATCTATATAGATATAGATATATATCTATATATAGATATAGATATATCTATATAGATATATATCTATAGATATCTATATATATAGATATATAGATATCTATATCTATATT",
            )

