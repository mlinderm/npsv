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
##contig=<ID=1,length=249250621>
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
##contig=<ID=2,length=243199373>
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


class SimpleDELVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	899922	HG3_Ill_SVrefine2DISCOVARDovetail_2	GGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGCTGTGGGAGGGGTCCGCGCGTCCGCAGTGGGGATGT	G	20	PASS	END=899992;SVTYPE=DEL;SVLEN=-70
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_event_length(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertEqual(variant.event_length, 70)

    def test_region_strings(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertEqual(variant.region_string(),"1:899923-899992")
        self.assertEqual(variant.left_flank_region_string(left_flank=2, right_flank=5),"1:899921-899927")
        self.assertEqual(variant.right_flank_region_string(left_flank=2, right_flank=5),"1:899991-899997")

    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertIsNotNone(variant)
        
        self.assertEqual(variant.left_flank_region_string(left_flank=1, right_flank=1), "1:899922-899923")
        self.assertEqual(variant.right_flank_region_string(left_flank=1, right_flank=1), "1:899992-899993")

        self.assertEqual(variant.ref_breakpoints(self.args.flank), ("1:1-2", "1:71-72"))
        self.assertEqual(variant.alt_breakpoints(self.args.flank), ("1:1-2", None))

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

        covg_path, ref_contig, alt_contig = variant.gnomad_coverage_profile(
            self.args, 
            os.path.join(FILE_DIR, "1_896922_903086.gnomad.genomes.coverage.summary.tsv.gz"),
            line_width=sys.maxsize,
        )
        self.assertEqual(ref_contig, "1_899922_899993")
        self.assertEqual(alt_contig, "1_899922_899993_alt")
        with open(covg_path, "r") as fasta:
            lines = [line.strip() for line in fasta]
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[1], "1_899922_899993_alt\t=1")


class SymbolicDELVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
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


class ComplexDELVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=4,length=191154276>
#CHROM POS ID REF ALT QUAL FILTER INFO
4	20473846	HG3_PB_assemblyticsfalcon_19066	TATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATA    ATATATATAGATATATCTATATATATCTATATAGATATATCTATATCTATATAGATATATCTATATATATATAGATATATCTATATCTATATAGATATATATCTATATATATATCTATATAGATATATCTATATAGATATAGATATATATCTATATATAGATATAGATATATCTATATAGATATATATCTATAGATATCTATATATATAGATATATAGATATCTATATCTATAT	10   LongReadHomRef	END=20474269;SVTYPE=DEL;SVLEN=-190
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def  test_variant_properties(self):
        record = next(vcf.Reader(self.vcf_file))
        self.assertTrue(record.is_sv)
        variant = Variant.from_pyvcf(record, None)
        self.assertIsNotNone(variant)
        self.assertTrue(variant.is_deletion)

        self.assertEqual(
            variant.ref_length, len("TATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATA")
        )

        self.assertEqual(
            variant.alt_length,
            len(
                "ATATATATAGATATATCTATATATATCTATATAGATATATCTATATCTATATAGATATATCTATATATATATAGATATATCTATATCTATATAGATATATATCTATATATATATCTATATAGATATATCTATATAGATATAGATATATATCTATATATAGATATAGATATATCTATATAGATATATATCTATAGATATCTATATATATAGATATATAGATATCTATATCTATAT"
            ),
        )

    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        
        self.assertEqual(variant.left_flank_region_string(left_flank=1, right_flank=1), "4:20473845-20473846")
        self.assertEqual(variant.right_flank_region_string(left_flank=1, right_flank=1), "4:20474269-20474270")

        self.assertEqual(variant.ref_breakpoints(self.args.flank), ("4:1-2", "4:425-426"))
        self.assertEqual(variant.alt_breakpoints(self.args.flank), ("4:1-2", "4:235-236"))

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GTATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATAT",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            mock_ref.assert_called_once_with(region="4:20473845-20474270")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">4_20473845_20474270")
            self.assertEqual(lines[0], f">{ref_contig}")
            self.assertEqual(
                lines[1],
                "GTATATATATATAGATCTATATATCTATATATAGATCTATATATAGATATATATCTATATATATAGATATATAGATATATAGATCTATATATAGATATATATATCTATATATAGATCTATATATAGATATAGATATCTATATAGATATCTATATCTATATATATGTAGATATATAGATATAGATATCTATATATCTATATATATAGATATCTATAGATATATATCTATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATATCTATATATAGATAGATATATATCTATATATAGATATATCTATATCTATATATAGATATATATCTATATATAGATATATCTATATATAGATATATATCTATAGATATATCTATATATATCGATATATCTATATATATCGATATATAT",
            )
            self.assertEqual(lines[2], ">4_20473845_20474270_alt")
            self.assertEqual(lines[2], f">{alt_contig}")
            self.assertEqual(
                lines[3],
                "GATATATATAGATATATCTATATATATCTATATAGATATATCTATATCTATATAGATATATCTATATATATATAGATATATCTATATCTATATAGATATATATCTATATATATATCTATATAGATATATCTATATAGATATAGATATATATCTATATATAGATATAGATATATCTATATAGATATATATCTATAGATATCTATATATATAGATATATAGATATCTATATCTATATT",
            )

class SingleBaseComplexDELVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=8,length=146364022>
#CHROM POS ID REF ALT QUAL FILTER INFO
8	79683398	.	AACCTCCCAACGCAATAGACATTGTGGTTTTCATTGCATATCATTCCTATTTCTCTCTCTCCATTATTTAGCAGTAATTTTTTTAATGAA	C	20	PASS	END=79683487;SVTYPE=DEL;SVLEN=-89
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        
        self.assertEqual(variant.left_flank_region_string(left_flank=1, right_flank=1), "8:79683397-79683398")
        self.assertEqual(variant.right_flank_region_string(left_flank=1, right_flank=1), "8:79683487-79683488")

        self.assertEqual(variant.ref_breakpoints(self.args.flank), ("8:1-2", "8:91-92"))
        self.assertEqual(variant.alt_breakpoints(self.args.flank), ("8:1-2", "8:2-3"))

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="AAACCTCCCAACGCAATAGACATTGTGGTTTTCATTGCATATCATTCCTATTTCTCTCTCTCCATTATTTAGCAGTAATTTTTTTAATGAAA"
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            mock_ref.assert_called_once_with(region="8:79683397-79683488")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">8_79683397_79683488")
            self.assertEqual(lines[0], f">{ref_contig}")
            self.assertEqual(
                lines[1],
                "AAACCTCCCAACGCAATAGACATTGTGGTTTTCATTGCATATCATTCCTATTTCTCTCTCTCCATTATTTAGCAGTAATTTTTTTAATGAAA",
            )
            self.assertEqual(lines[2], ">8_79683397_79683488_alt")
            self.assertEqual(lines[2], f">{alt_contig}")
            self.assertEqual(
                lines[3],
                "ACA",
            )

class NonNormalizedComplexDELVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=4,length=191154276>
#CHROM POS ID REF ALT QUAL FILTER INFO
4	32197281	.	TGAACCTGGGAGGCAGAGCTTGCAGTGAGCAGAGATCATGCCACTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAATTAGCCAGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATGGCATG	TG	.	PASS    END=32197447;SVTYPE=DEL;SVLEN=-165
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()    
    
    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        
        self.assertEqual(variant.left_flank_region_string(left_flank=1, right_flank=1), "4:32197282-32197283")
        self.assertEqual(variant.right_flank_region_string(left_flank=1, right_flank=1), "4:32197447-32197448")

        self.assertEqual(variant.ref_breakpoints(self.args.flank), ("4:1-2", "4:166-167"))
        self.assertEqual(variant.alt_breakpoints(self.args.flank), ("4:1-2", None))

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="GAACCTGGGAGGCAGAGCTTGCAGTGAGCAGAGATCATGCCACTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAATTAGCCAGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATGGCATGA"
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            mock_ref.assert_called_once_with(region="4:32197282-32197448")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">4_32197282_32197448")
            self.assertEqual(lines[0], f">{ref_contig}")
            self.assertEqual(
                lines[1],
                "GAACCTGGGAGGCAGAGCTTGCAGTGAGCAGAGATCATGCCACTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAATTAGCCAGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATGGCATGA",
            )
            self.assertEqual(lines[2], ">4_32197282_32197448_alt")
            self.assertEqual(lines[2], f">{alt_contig}")
            self.assertEqual(
                lines[3],
                "GA",
            )

class SimpleINSVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=INS,Description="Insertion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	931634	HG2_PB_SVrefine2PB10Xhap12_17	A	AGGGAGGGCAGAAAGGACCCCCACGTGAGGGGGCACCCCACATCTGGGGCCACAGGATGCAGGGTGGGGAGGGCAGAAAGGCCCCCCCGCGGGAAGGGGCACCCCACATCTGGGCCACAGGATGCAGGGTGGGGAGGGCAGAAAGGCCCCCCCGCGGGAAGGGGCACCCCACATCTGGGGCCACAGGATGCAGGGTG	.	PASS	SVTYPE=INS;END=931634;SVLEN=196
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_region_string(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        with self.assertRaises(ValueError):
            variant.region_string()
        self.assertEqual(variant.region_string(flank=1), "1:931634-931635")

    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertIsNotNone(variant)
        
        self.assertEqual(variant.left_flank_region_string(left_flank=1, right_flank=1), "1:931634-931635")
        self.assertEqual(variant.right_flank_region_string(left_flank=1, right_flank=1), "1:931634-931635")

        self.assertEqual(variant.ref_breakpoints(self.args.flank), ("1:1-2", None))
        self.assertEqual(variant.alt_breakpoints(self.args.flank), ("1:1-2", "1:197-198"))

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="AG",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            self.assertEqual(ref_contig, "1_931634_931635")
            self.assertEqual(alt_contig, "1_931634_931635_alt")
            mock_ref.assert_called_once_with(region="1:931634-931635")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">1_931634_931635")
            self.assertEqual(lines[1], "AG")
            self.assertEqual(lines[2], ">1_931634_931635_alt")
            self.assertEqual(lines[3], "AGGGAGGGCAGAAAGGACCCCCACGTGAGGGGGCACCCCACATCTGGGGCCACAGGATGCAGGGTGGGGAGGGCAGAAAGGCCCCCCCGCGGGAAGGGGCACCCCACATCTGGGCCACAGGATGCAGGGTGGGGAGGGCAGAAAGGCCCCCCCGCGGGAAGGGGCACCCCACATCTGGGGCCACAGGATGCAGGGTGG")

class SymbolicINSVariantTestSuite(unittest.TestCase):
    def test_manta_vcf(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description="Length of insertion">
##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion">
##ALT=<ID=INS,Description="Insertion">
##contig=<ID=3,length=198022430>
#CHROM POS ID REF ALT QUAL FILTER INFO
3       72386664        MantaINS:1:5470:5470:0:0:0      G       <INS>   999     MaxDepth        END=72386664;SVTYPE=INS;SVLEN=10;CIPOS=0,9;CIEND=0,9;SVINSLEN=10;SVINSSEQ=GTGTGTGTGC
"""
        )
        record = next(vcf.Reader(vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertEqual(variant._alt_seq(flank=1,ref_seq="GT"), "GGTGTGTGTGCT")

class SymbolicDUPVariantTestSuite(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DUP,Description="Duplication">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1   4999478 19  N   <DUP>   4133.49 PASS    SVTYPE=DUP;SVLEN=254;END=4999732
"""
        )
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)
    
    def tearDown(self):
        self.tempdir.cleanup()

    def test_region_string(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertIsNotNone(variant)
        self.assertEqual(variant.region_string(flank=1), "1:4999478-4999733")

    def test_breakpoints(self):
        record = next(vcf.Reader(self.vcf_file))
        variant = Variant.from_pyvcf(record, None)
        self.assertIsNotNone(variant)
        
        self.assertEqual(variant.left_flank_region_string(left_flank=1, right_flank=1), "1:4999478-4999479")
        self.assertEqual(variant.right_flank_region_string(left_flank=1, right_flank=1), "1:4999732-4999733")

        self.assertEqual(variant.ref_breakpoints(self.args.flank), ("1:1-2", "1:255-256"))
        self.assertEqual(variant.alt_breakpoints(self.args.flank), ("1:255-256", None))
        

    def test_consensus_fasta(self):
        with patch.object(
            Variant,
            "reference_sequence",
            return_value="TCTCCATATGATGTCAGTGTCCTCCATATGATGTCAGTGTCCTCCATATGACATCAATATCCTCCATATGATATCAATATCCTCTGTATTGATATTGATATTGATATTTGGAGGATATCAATATCCTCCAAATGATGTCAGTGTCCTCCATATGATGTCAATGTCCTCCATATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCATATGATGTCAGTGTCCTCTGTATGAC",
        ) as mock_ref:
            record = next(vcf.Reader(self.vcf_file))
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

            fasta_path, ref_contig, alt_contig = variant.synth_fasta(self.args, line_width=sys.maxsize)
            mock_ref.assert_called_once_with(region="1:4999478-4999733")
            
            self.assertEqual(ref_contig, "1_4999478_4999733")
            self.assertEqual(alt_contig, "1_4999478_4999733_alt")

            with open(fasta_path, "r") as fasta:
                lines = [line.strip() for line in fasta]
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], ">1_4999478_4999733")
            self.assertEqual(lines[1], "TCTCCATATGATGTCAGTGTCCTCCATATGATGTCAGTGTCCTCCATATGACATCAATATCCTCCATATGATATCAATATCCTCTGTATTGATATTGATATTGATATTTGGAGGATATCAATATCCTCCAAATGATGTCAGTGTCCTCCATATGATGTCAATGTCCTCCATATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCATATGATGTCAGTGTCCTCTGTATGAC")
            self.assertEqual(lines[2], ">1_4999478_4999733_alt")
            self.assertEqual(lines[3], "TCTCCATATGATGTCAGTGTCCTCCATATGATGTCAGTGTCCTCCATATGACATCAATATCCTCCATATGATATCAATATCCTCTGTATTGATATTGATATTGATATTTGGAGGATATCAATATCCTCCAAATGATGTCAGTGTCCTCCATATGATGTCAATGTCCTCCATATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCATATGATGTCAGTGTCCTCTGTATGACTCCATATGATGTCAGTGTCCTCCATATGATGTCAGTGTCCTCCATATGACATCAATATCCTCCATATGATATCAATATCCTCTGTATTGATATTGATATTGATATTTGGAGGATATCAATATCCTCCAAATGATGTCAGTGTCCTCCATATGATGTCAATGTCCTCCATATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCGTATGATGTCAATATCCTCCATATGATGTCAGTGTCCTCTGTATGAC")
