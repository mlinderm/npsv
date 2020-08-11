import argparse, io, os, sys, tempfile, unittest
from unittest.mock import patch
import vcf
import pysam
from npsv.variant import Variant
from npsv.sample import Sample
from npsv.feature_extraction import extract_deletion_features
from npsv import npsva

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class DeletionFeaturesTest(unittest.TestCase):
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
        self.args = argparse.Namespace(flank=3000, min_anchor=11, default_ci=10, min_mapq=40, min_baseq=15, rel_coverage_flank=1000)
        self.input_fasta = os.path.join(FILE_DIR, "1_2073761_2073846_DEL.synth.fasta")
        self.input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
        self.sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), self.input_bam)

        patcher = patch.object(Variant, "reference_sequence", return_value="AGCAGCCGAAGCGCCTCCTTTCAATCCAGGGTCCACACATCCAGCAGCCGAAGCGCCCTCCTTTCAATCCAGGGTCCAGGCATCT")
        self.mock_reference_sequence = patcher.start()
        self.addCleanup(patcher.stop)

    def test_read_gather_call(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            
            with patch.object(npsva.RealignedFragments, "gather_reads", return_value=0) as mock_gather:
                variant = Variant.from_pyvcf(record, None)
                extract_deletion_features(
                    self.args,
                    variant,
                    self.input_bam,
                    self.sample,
                    input_fasta=self.input_fasta,
                    ref_contig="1_2073761_2073846_DEL",
                    alt_contig="1_2073761_2073846_DEL_alt",
                )
                # Since this is a small variant, there should only one gather reads call
                mock_gather.assert_called_once()
    
    def test_feature_extraction(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
    
            variant = Variant.from_pyvcf(record, None)
            features = extract_deletion_features(
                self.args,
                variant,
                self.input_bam,
                self.sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )
            features.print_features(sys.stdout)
    