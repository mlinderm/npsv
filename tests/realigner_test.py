import argparse, io, os, tempfile, unittest
from unittest.mock import patch
import vcf
import pysam
from npsv.variant import Variant
from npsv.feature_extraction import count_realigned_reads
from npsv.sample import Sample
from npsv import npsva

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class DELBreakpoints(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)
        self.input_bam = "dummy.bam"
        self.sample = Sample.from_npsv(
            os.path.join(FILE_DIR, "stats.json"), self.input_bam
        )

    def tearDown(self):
        self.tempdir.cleanup()

    @patch('npsv.npsva.RealignedFragments')
    def test_resolved_breakpoint_arguments(self, mock_fragments):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	100	.	ACGT	A	.	PASS	SVTYPE=DEL;END=103;SVLEN=-3
"""
        )

        for record in vcf.Reader(vcf_file):
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

        read_counts = {"rl": 0, "rr": 0, "al": 0, "ar": 0}
        mock_fragments.count_realigned_reads.return_value = (read_counts, [])

        count_realigned_reads(self.args, mock_fragments, variant)
        mock_fragments.count_realigned_reads.assert_called_once_with(
            [("ref:1-2","ref:4-5","alt:1-2","")],
            count_straddle=True,
        )

    @patch('npsv.npsva.RealignedFragments')
    def test_symbolic_breakpoint_arguments(self, mock_fragments):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	100	.	A	<DEL>	.	PASS	SVTYPE=DEL;END=105;SVLEN=-5
"""
        )

        for record in vcf.Reader(vcf_file):
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

        read_counts = {"rl": 0, "rr": 0, "al": 0, "ar": 0}
        mock_fragments.count_realigned_reads.return_value = (read_counts, [])

        count_realigned_reads(self.args, mock_fragments, variant)
        mock_fragments.count_realigned_reads.assert_called_once_with(
            [("ref:1-2","ref:6-7","alt:1-2","")],
            count_straddle=True,
        )

    @patch('npsv.npsva.RealignedFragments')
    def test_complex_breakpoint_arguments(self, mock_fragments):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	100	.	AACGGTT	ACGT	.	PASS	SVTYPE=DEL;END=106;SVLEN=-3
"""
        )

        for record in vcf.Reader(vcf_file):
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

        read_counts = {"rl": 0, "rr": 0, "al": 0, "ar": 0}
        mock_fragments.count_realigned_reads.return_value = (read_counts, [])

        count_realigned_reads(self.args, mock_fragments, variant)
        mock_fragments.count_realigned_reads.assert_called_once_with(
            [("ref:1-2","ref:7-8","alt:1-2","alt:4-5")],
            count_straddle=True,
        )


class INSBreakpoints(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)
        self.input_bam = "dummy.bam"
        self.sample = Sample.from_npsv(
            os.path.join(FILE_DIR, "stats.json"), self.input_bam
        )

    def tearDown(self):
        self.tempdir.cleanup()

    @patch('npsv.npsva.RealignedFragments')
    def test_resolved_breakpoint_arguments(self, mock_fragments):
        self.vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	100	.	A	ACGT	.	PASS	SVTYPE=INS;END=100;SVLEN=3
"""
        )
        for record in vcf.Reader(self.vcf_file):
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

        read_counts = {"rl": 0, "rr": 0, "al": 0, "ar": 0}
        mock_fragments.count_realigned_reads.return_value = (read_counts, [])

        count_realigned_reads(self.args, mock_fragments, variant)
        mock_fragments.count_realigned_reads.assert_called_once_with(
            [("ref:1-2","","alt:1-2","alt:4-5")],
            count_straddle=True,
        )


class NPSVARealignedFragmentsTest(unittest.TestCase):
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

    def test_straddle_implementation(self):
        try:
            # Create SAM file with a single read
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".sam"
            ) as sam_file:
                # fmt: off
                print("@HD", "VN:1.3", "SO:coordinate", sep="\t", file=sam_file)
                print("@SQ", "SN:ref", "LN:6070", sep="\t", file=sam_file)
                print("@RG", "ID:synth1", "LB:synth1", "PL:illumina", "PU:ART", "SM:HG002", sep="\t", file=sam_file)
                print(
                    "ref-82", "99", "ref", "91", "99", "148=", "=", "679", "736",
                    "GATGAGCGAGAGCCGCCAGACCCACGTGACGCTGCACGACATCGACCCTCAGGCCTTGGACCAGCTGGTGCAGTTTGCCTACACGGCTGAGATTGTGGTGGGCGAGGGCAATGTGCAGGTGAGGGCTCCCTCACCCGGATCCCGGTGT",
                    "CCCGGGGGGGGGG=CGJJGCJJJJJJJJJGJJJGCGGJJJJJJGJJGJCG8GGJGJJJGGCGCJGCCJCCGGG81GGCGGGGCCCGGCGGGGGGGC=GCCGGCGGCCGGGCCGGGC8CGGGCCC=GGCGGGGGGGGGGGCGGGGGGCG",
                    sep="\t", file=sam_file,
                )
                print(
                    "ref-82", "147", "ref", "679", "99", "148=", "=", "91", "-736",
                    "CCTGACTCTGCTCGGCCCCTCCCAGTATGAACACTCAGCCCCCACCTGCTAACCCTCCCTCCTAGGCATCTTCAGGGCTCCCTGGGTCCACAGGACCCTCCCCAGATCTCAGGTCTGAGGACCCCCACTCCCAGGTTCTGGAACTGGT",
                    "CCGGCGGGGCCCGCC=GGGG8CG8GGGC8CCGGGGGGGGGGGGGCCC(JGG1CGGGGCGGCCGC8GGGCGGGGGCCGGGJCGGG(CJ=JGJJGJJGGJJGGJJGGGJGJJGJJGJJGJGCCCGJGJJGJJJJJGJGGCG1GGGGGCCC",
                    sep="\t", file=sam_file,
                )
                # fmt: on

            # Strict straddler
            self.assertTrue(npsva.test_straddle(sam_file.name, "ref:1-101", "ref:588-688", 10, True))
            self.assertFalse(npsva.test_straddle(sam_file.name, "ref:1-101", "ref:817-917", 10, True))

            # Not-strict straddler
            self.assertTrue(npsva.test_straddle(sam_file.name, "ref:1-101", "ref:817-917", 10, False))

        finally:
            os.remove(sam_file.name)

    def test_pipeline_straddle_counting(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), input_bam)

            fragments = npsva.RealignedFragments(
                self.input_fasta,
                sample.mean_insert_size,
                sample.std_insert_size,
                sample.insert_size_density().as_dict(),
                input_bam,
            )
            fragments.gather_reads(variant.region_string(flank=self.args.flank))
            self.assertEqual(fragments.size(), 254)

            left_breakpoint = variant.left_flank_region_string(left_flank=1, right_flank=1)
            right_breakpoint = variant.right_flank_region_string(left_flank=1, right_flank=1)
            pair_results = fragments.count_pipeline_straddlers(
                left_breakpoint, right_breakpoint, self.args.flank, -variant.event_length, 1.5, 10,
            )
            self.assertAlmostEqual(pair_results["alt_weighted_count"], 13.496, places=1)
            self.assertAlmostEqual(pair_results["insert_lower"], 0.0, places=2)
            self.assertAlmostEqual(pair_results["insert_upper"] / pair_results["insert_count"], 0.166, places=2)

    def test_pipeline_clip_counting(self):
        try:
            # Create SAM file with a single read
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".sam"
            ) as sam_file:
                # fmt: off
                print("@HD", "VN:1.3", "SO:coordinate", sep="\t", file=sam_file)
                print("@SQ", "SN:1", "LN:6070", sep="\t", file=sam_file)
                print("@RG", "ID:synth1", "LB:synth1", "PL:illumina", "PU:ART", "SM:HG002", sep="\t", file=sam_file)
                print(
                    "ref-82", "99", "1", "91", "99", "148=", "=", "679", "732",
                    "GATGAGCGAGAGCCGCCAGACCCACGTGACGCTGCACGACATCGACCCTCAGGCCTTGGACCAGCTGGTGCAGTTTGCCTACACGGCTGAGATTGTGGTGGGCGAGGGCAATGTGCAGGTGAGGGCTCCCTCACCCGGATCCCGGTGT",
                    "CCCGGGGGGGGGG=CGJJGCJJJJJJJJJGJJJGCGGJJJJJJGJJGJCG8GGJGJJJGGCGCJGCCJCCGGG81GGCGGGGCCCGGCGGGGGGGC=GCCGGCGGCCGGGCCGGGC8CGGGCCC=GGCGGGGGGGGGGGCGGGGGGCG",
                    sep="\t", file=sam_file,
                )
                print(
                    "ref-82", "147", "1", "679", "99", "4S140=4H", "=", "91", "-732",
                    "CCTGACTCTGCTCGGCCCCTCCCAGTATGAACACTCAGCCCCCACCTGCTAACCCTCCCTCCTAGGCATCTTCAGGGCTCCCTGGGTCCACAGGACCCTCCCCAGATCTCAGGTCTGAGGACCCCCACTCCCAGGTTCTGGAAC",
                    "CCGGCGGGGCCCGCC=GGGG8CG8GGGC8CCGGGGGGGGGGGGGCCC(JGG1CGGGGCGGCCGC8GGGCGGGGGCCGGGJCGGG(CJ=JGJJGJJGGJJGGJJGGGJGJJGJJGJJGJGCCCGJGJJGJJJJJGJGGCG1GGGG",
                    sep="\t", file=sam_file,
                )
                # fmt: on

            sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), sam_file.name)
            fragments = npsva.RealignedFragments(
                self.input_fasta,
                sample.mean_insert_size,
                sample.std_insert_size,
                sample.insert_size_density().as_dict(),
                sam_file.name,
            )
            fragments.gather_reads("1:1-1000")
            self.assertEqual(fragments.size(), 1)
            
            clip_results = fragments.count_pipeline_clipped_reads("1:100-101", 4)
            self.assertDictEqual(clip_results, {"left": 0, "right": 0, "both": 0, "total": 1})

            clip_results = fragments.count_pipeline_clipped_reads("1:700-701", 4)
            self.assertDictEqual(clip_results, {"left": 0, "right": 0, "both": 1, "total": 1})

        finally:
            os.remove(sam_file.name)

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

    def test_breakpoint_overlap(self):
        try:
            # Create SAM file with a single read
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".sam"
            ) as sam_file:
                # fmt: off
                print("@HD", "VN:1.3", "SO:coordinate", sep="\t", file=sam_file)
                print("@SQ", "SN:ref", "LN:6070", sep="\t", file=sam_file)
                print("@RG", "ID:synth1", "LB:synth1", "PL:illumina", "PU:ART", "SM:HG002", sep="\t", file=sam_file)
                print(
                    "ref-82", "99", "ref", "91", "99", "148=", "=", "679", "736",
                    "GATGAGCGAGAGCCGCCAGACCCACGTGACGCTGCACGACATCGACCCTCAGGCCTTGGACCAGCTGGTGCAGTTTGCCTACACGGCTGAGATTGTGGTGGGCGAGGGCAATGTGCAGGTGAGGGCTCCCTCACCCGGATCCCGGTGT",
                    "CCCGGGGGGGGGG=CGJJGCJJJJJJJJJGJJJGCGGJJJJJJGJJGJCG8GGJGJJJGGCGCJGCCJCCGGG81GGCGGGGCCCGGCGGGGGGGC=GCCGGCGGCCGGGCCGGGC8CGGGCCC=GGCGGGGGGGGGGGCGGGGGGCG",
                    sep="\t", file=sam_file,
                )
                print(
                    "ref-82", "147", "ref", "679", "99", "148=", "=", "91", "-736",
                    "CCTGACTCTGCTCGGCCCCTCCCAGTATGAACACTCAGCCCCCACCTGCTAACCCTCCCTCCTAGGCATCTTCAGGGCTCCCTGGGTCCACAGGACCCTCCCCAGATCTCAGGTCTGAGGACCCCCACTCCCAGGTTCTGGAACTGGT",
                    "CCGGCGGGGCCCGCC=GGGG8CG8GGGC8CCGGGGGGGGGGGGGCCC(JGG1CGGGGCGGCCGC8GGGCGGGGGCCGGGJCGGG(CJ=JGJJGJJGGJJGGJJGGGJGJJGJJGJJGJGCCCGJGJJGJJJJJGJGGCG1GGGGGCCC",
                    sep="\t", file=sam_file,
                )
                # fmt: on

            # Doesn't overlap both bases, i.e. doesn't span breakpoint
            self.assertFalse(
                npsva.test_alignment_overlap(sam_file.name, "ref:90-91", False)
            )
            self.assertFalse(
                npsva.test_alignment_overlap(sam_file.name, "ref:826-827", False)
            )

            # Does overlap both bases, i.e. does span breakpoint
            self.assertTrue(
                npsva.test_alignment_overlap(sam_file.name, "ref:91-92", False)
            )
            self.assertTrue(
                npsva.test_alignment_overlap(sam_file.name, "ref:825-826", False)
            )

            # Positions in the insert between reads (with and without counting straddlers)
            self.assertTrue(
                npsva.test_alignment_overlap(sam_file.name, "ref:247-248", True)
            )
            self.assertFalse(
                npsva.test_alignment_overlap(sam_file.name, "ref:247-248", False)
            )

        finally:
            os.remove(sam_file.name)

    def test_realigned_read_counting(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), input_bam)

            fragments = npsva.RealignedFragments(
                self.input_fasta,
                sample.mean_insert_size,
                sample.std_insert_size,
                sample.insert_size_density().as_dict(),
                input_bam,
            )
            fragments.gather_reads(variant.region_string(flank=self.args.flank))
            self.assertEqual(fragments.size(), 254)

            ref_contig = "1_2073761_2073846_DEL"
            alt_contig = "1_2073761_2073846_DEL_alt"
            
            rl_breakpoint = f"{ref_contig}:{self.args.flank}-{self.args.flank+1}"
            al_breakpoint = f"{alt_contig}:{self.args.flank}-{self.args.flank+1}"
            ref_length = variant.ref_length    
            rr_breakpoint = f"{ref_contig}:{self.args.flank + ref_length - 1}-{self.args.flank + ref_length}"

            counts, read_names = fragments.count_realigned_reads([(rl_breakpoint, rr_breakpoint, al_breakpoint, "")])            
            self.assertEqual(counts["al"], 18.0)
            self.assertEqual((counts["rl"] + counts["rr"]) / 2, 4.0)
            for bp in ("rl", "rr", "al", "rl"):
                self.assertEqual(len(read_names[bp]), counts[bp])

    