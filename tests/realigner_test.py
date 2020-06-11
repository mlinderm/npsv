import argparse, io, os, tempfile, unittest
import vcf
import pysam
from npsv.feature_extraction import Variant
from npsv.sample import Sample
from npsv import npsva

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class NPSVAAlleleCountingTest(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1 2073761 . CAGCAGCCGAAGCGCCTCCTTTCAATCCAGGGTCCACACATCCAGCAGCCGAAGCGCCCTCCTTTCAATCCAGGGTCCAGGCATCT C . PASS SVTYPE=DEL;END=2073846;SVLEN=-85
"""
        )
        self.args = argparse.Namespace(flank=3000)
        self.input_fasta = os.path.join(FILE_DIR, "1_2073761_2073846_DEL.synth.fasta")

    def test_count_alleles(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_distribution(
                input_bam, 569.2, 163, 148, mean_coverage=25
            )
            hom_alt_ref, hom_alt_alt = variant.count_alleles_with_npsva(
                self.args,
                input_bam,
                sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )
            self.assertEqual(hom_alt_ref, 3.0)
            self.assertEqual(hom_alt_alt, 16.0)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_1.bam")
            sample = Sample.from_distribution(
                input_bam, 569.2, 163, 148, mean_coverage=25
            )
            het_ref, het_alt = variant.count_alleles_with_npsva(
                self.args,
                input_bam,
                sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )
            self.assertEqual(het_ref, 16.0)
            self.assertEqual(het_alt, 21.0)

    def test_overlap_bam(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_distribution(
                input_bam, 569.2, 163, 148, mean_coverage=25
            )

            try:
                output_bam_file = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".bam"
                )
                output_bam_file.close()
                hom_alt_ref, hom_alt_alt = variant.count_alleles_with_npsva(
                    self.args,
                    input_bam,
                    sample,
                    input_fasta=self.input_fasta,
                    ref_contig="1_2073761_2073846_DEL",
                    alt_contig="1_2073761_2073846_DEL_alt",
                    overlap_bam=output_bam_file.name,
                )
                self.assertTrue(os.path.exists(output_bam_file.name))

                bam_reader = pysam.AlignmentFile(output_bam_file.name)
                read_count = 0
                for read in bam_reader.fetch():
                    read_count += 1
                    self.assertTrue(read.has_tag("ov"))
                    self.assertTrue(read.has_tag("as"))
                    self.assertTrue(read.has_tag("fs"))
                self.assertGreaterEqual(
                    read_count, hom_alt_alt
                )  # At least one read for each alternate fragment

            finally:
                os.remove(output_bam_file.name)

    def test_alignment_scoring(self):
        try:
            # Create SAM file with a single read
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".sam"
            ) as sam_file:
                # fmt: off
                print("@HD", "VN:1.3", "SO:coordinate", sep="\t", file=sam_file)
                print("@SQ", "SN:1", "LN:249250621", sep="\t", file=sam_file)
                print("@RG", "ID:synth1", "LB:synth1", "PL:illumina", "PU:ART", "SM:HG002", sep="\t", file=sam_file)
                print(
                    "ref-354", "147", "1", "1", "60", "148M", "=", "2073433", "-435",
                    "AGCAGCCGAAGCGCCTCCTTTCACTCTAGGGTCCAGGCATCCAGCAGCCGAAGCGCCTCCTTTCAATCCAGGGTCCACACATCCAGCAGCCGAAGCGCCCTCCTTTCAATCCAGGGTCCAGGCATCTAGCAGCCGAAGCGCCTCCTTT",
                    "GG8CCGGGCGGGCGGGGGCGGCGGGGGGGGGGGCGGCGGGG=GGGJCCJGGGGGGCGGGGGGCG1GGCGG8GGCGC1GGCGJGCCGGJGJGJGGCGCJGJGJJCGGJJCJJGJJGJJJGJGCJJJGGJJJJGJJJGGGCGGGCGGCCC",
                    "RG:Z:synth1",
                    sep="\t",
                    file=sam_file,
                )
                # fmt: on

            # Read was aligned to very beginning of reference, so using read as reference should be all matches
            ref_sequence = "AGCAGCCGAAGCGCCTCCTTTCACTCTAGGGTCCAGGCATCCAGCAGCCGAAGCGCCTCCTTTCAATCCAGGGTCCACACATCCAGCAGCCGAAGCGCCCTCCTTTCAATCCAGGGTCCAGGCATCTAGCAGCCGAAGCGCCTCCTTT"
            scores = npsva.test_score_alignment(ref_sequence, sam_file.name)
            self.assertEqual(len(scores), 1)
            self.assertLess(scores[0], 0)
        finally:
            os.remove(sam_file.name)
