import argparse, io, os, unittest
import vcf
from npsv.feature_extraction import Variant
from npsv.sample import Sample

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class NPSVAAlleleCountingTest(unittest.TestCase):
    def test_count_alleles(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1 2073761 . CAGCAGCCGAAGCGCCTCCTTTCAATCCAGGGTCCACACATCCAGCAGCCGAAGCGCCCTCCTTTCAATCCAGGGTCCAGGCATCT C . PASS SVTYPE=DEL;END=2073846;SVLEN=-85
"""
        )

        args = argparse.Namespace(flank=3000)
        input_fasta = os.path.join(FILE_DIR, "1_2073761_2073846_DEL.synth.fasta")

        for record in vcf.Reader(vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_distribution(input_bam, 569.2, 163, 148, mean_coverage=25)
            hom_alt_ref, hom_alt_alt = variant.count_alleles_with_npsva(
                args,
                input_bam,
                sample,
                input_fasta=input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_1.bam")
            sample = Sample.from_distribution(input_bam, 569.2, 163, 148, mean_coverage=25)
            het_ref, het_alt = variant.count_alleles_with_npsva(
                args,
                input_bam,
                sample,
                input_fasta=input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )
            self.assertGreater(hom_alt_alt / (hom_alt_ref + hom_alt_alt), het_alt / (het_ref + het_alt))
