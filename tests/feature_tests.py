import argparse, io, os, sys, tempfile, unittest
from unittest.mock import patch
import vcf
import pysam
from npsv.variant import Variant
from npsv.sample import Sample
from npsv.feature_extraction import extract_features, count_realigned_reads
from npsv import npsva

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class DELFeaturesTest(unittest.TestCase):
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
        self.args = argparse.Namespace(flank=3000, min_anchor=11, default_ci=10, min_mapq=40, min_baseq=15, rel_coverage_flank=1000, count_straddle=True, min_clip=4, mapq_reads=False, reference=None)
        self.input_fasta = os.path.join(FILE_DIR, "1_2073761_2073846_DEL.synth.fasta")
        self.input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
        self.sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), self.input_bam)

        patcher = patch.object(Variant, "reference_sequence", return_value="AGCAGCCGAAGCGCCTCCTTTCAATCCAGGGTCCACACATCCAGCAGCCGAAGCGCCCTCCTTTCAATCCAGGGTCCAGGCATCT")
        self.mock_reference_sequence = patcher.start()
        self.addCleanup(patcher.stop)

    def test_read_gather_call(self):
        for record in vcf.Reader(self.vcf_file):
            with patch.object(npsva.RealignedFragments, "gather_reads", return_value=0) as mock_gather:
                variant = Variant.from_pyvcf(record, None)
                extract_features(
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
    
    def test_count_realigned_alleles(self):
        for record in vcf.Reader(self.vcf_file):
            variant = Variant.from_pyvcf(record, None)
            
            fragments = npsva.RealignedFragments(
                self.input_fasta,
                self.sample.mean_insert_size,
                self.sample.std_insert_size,
                self.sample.insert_size_density().as_dict(),
                self.input_bam,
            )
            fragments.gather_reads(variant.region_string(flank=self.args.flank))
            self.assertGreater(fragments.size(), 0)

            hom_alt_ref, hom_alt_alt, _ = count_realigned_reads(
                self.args,
                fragments,
                variant,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )
            self.assertTrue(hom_alt_ref, 4.0)
            self.assertTrue(hom_alt_alt, 18.0)


    def test_feature_extraction(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
    
            variant = Variant.from_pyvcf(record, None)
            features = extract_features(
                self.args,
                variant,
                self.input_bam,
                self.sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )

@unittest.skipIf(not os.path.exists("/data/human_g1k_v37.fasta"), "Reference genome not available")
class DELSpanningFeaturesTest(unittest.TestCase):
    def setUp(self):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM POS ID REF ALT QUAL FILTER INFO
1	899922	.	G	<DEL>	.	.	SVTYPE=DEL;END=899998;SVLEN=-76;ORIGINAL=HG3_Ill_SVrefine2DISCOVARDovetail_2
"""
        )
        self.args = argparse.Namespace(flank=3000, min_anchor=11, default_ci=10, min_mapq=40, min_baseq=15, rel_coverage_flank=1000, count_straddle=True, min_clip=4, mapq_reads=False, reference="/data/human_g1k_v37.fasta", tempdir=tempfile.gettempdir())
        self.input_bam = os.path.join(FILE_DIR, "1_896922_902998.bam")
        self.sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), self.input_bam)

    def test_spanning_features(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
    
            variant = Variant.from_pyvcf(record, self.args.reference)
            features = extract_features(
                self.args,
                variant,
                self.input_bam,
                self.sample,
            )
            

class PySAMFeatures(unittest.TestCase):
    def test_region_parsing(self):
        # parse_region converts region to 0-indexed half-open coordinates
        self.assertEqual(pysam.libcutils.parse_region(region="1:1-1"), ('1', 0, 1))