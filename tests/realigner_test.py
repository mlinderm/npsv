import argparse, io, os, tempfile, unittest
from unittest.mock import patch
import vcf
import pysam
from npsv.variant import Variant
from npsv.feature_extraction import count_alleles_with_npsva
from npsv.sample import Sample
from npsv.fragment import SpanningFragments, gather_reads
from npsv import npsva

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")


class DELAlleleCountingTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)
        self.input_bam = "dummy.bam"
        self.sample = Sample.from_npsv(
            os.path.join(FILE_DIR, "stats.json"), self.input_bam
        )

    def tearDown(self):
        self.tempdir.cleanup()

    def test_resolved_breakpoint_arguments(self):
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

        read_counts = {"rl_reads": 0, "rr_reads": 0, "al_reads": 0, "ar_reads": 0}

        with patch.object(
            Variant, "reference_sequence", return_value="AT"
        ), patch.object(
            npsva.Realigner, "count_alignments", return_value=(read_counts, [])
        ) as mock_count:
            count_alleles_with_npsva(
                self.args, variant, self.input_bam, self.sample,
            )
            mock_count.assert_called_once_with(
                self.input_bam,
                "1_100_104:1-2",
                "1_100_104_alt:1-2",
                count_straddle=True,
                region="1:100-104",
                rr_region="1_100_104:4-5",
            )

    def test_symbolic_breakpoint_arguments(self):
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

        read_counts = {"rl_reads": 0, "rr_reads": 0, "al_reads": 0, "ar_reads": 0}

        with patch.object(
            Variant, "reference_sequence", return_value="AT"
        ), patch.object(
            npsva.Realigner, "count_alignments", return_value=(read_counts, [])
        ) as mock_count:
            count_alleles_with_npsva(
                self.args, variant, self.input_bam, self.sample,
            )
            mock_count.assert_called_once_with(
                self.input_bam,
                "1_100_106:1-2",
                "1_100_106_alt:1-2",
                count_straddle=True,
                region="1:100-106",
                rr_region="1_100_106:6-7",
            )

    def test_complex_breakpoint_arguments(self):
        vcf_file = io.StringIO(
            """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
#CHROM POS ID REF ALT QUAL FILTER INFO
1	100	.	ACCGGTT	ACGT	.	PASS	SVTYPE=DEL;END=106;SVLEN=-3
"""
        )

        for record in vcf.Reader(vcf_file):
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

        read_counts = {"rl_reads": 0, "rr_reads": 0, "al_reads": 0, "ar_reads": 0}

        with patch.object(
            Variant, "reference_sequence", return_value="AT"
        ), patch.object(
            npsva.Realigner, "count_alignments", return_value=(read_counts, [])
        ) as mock_count:
            count_alleles_with_npsva(
                self.args, variant, self.input_bam, self.sample,
            )
            mock_count.assert_called_once_with(
                self.input_bam,
                "1_100_107:1-2",
                "1_100_107_alt:1-2",
                ar_region="1_100_107_alt:4-5",
                count_straddle=True,
                region="1:100-107",
                rr_region="1_100_107:7-8",
            )


class INSAlleleCountingTest(unittest.TestCase):
    def setUp(self):
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
        self.tempdir = tempfile.TemporaryDirectory()
        self.args = argparse.Namespace(flank=1, tempdir=self.tempdir.name)
        self.input_bam = "dummy.bam"
        self.sample = Sample.from_npsv(
            os.path.join(FILE_DIR, "stats.json"), self.input_bam
        )

    def tearDown(self):
        self.tempdir.cleanup()

    def test_breakpoint_arguments(self):
        for record in vcf.Reader(self.vcf_file):
            variant = Variant.from_pyvcf(record, None)
            self.assertIsNotNone(variant)

        read_counts = {"rl_reads": 0, "rr_reads": 0, "al_reads": 0, "ar_reads": 0}

        with patch.object(
            Variant, "reference_sequence", return_value="AT"
        ), patch.object(
            npsva.Realigner, "count_alignments", return_value=(read_counts, [])
        ) as mock_count:
            count_alleles_with_npsva(
                self.args, variant, self.input_bam, self.sample,
            )
            mock_count.assert_called_once_with(
                self.input_bam,
                "1_100_101:1-2",
                "1_100_101_alt:1-2",
                ar_region="1_100_101_alt:4-5",
                count_straddle=True,
                region="1:100-101",
            )


class OverlapCoordinateTest(unittest.TestCase):
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
            variant = Variant.from_pyvcf(record, None)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_distribution(
                input_bam, 569.2, 163, 148, mean_coverage=25
            )
            hom_alt_ref, hom_alt_alt, read_names = count_alleles_with_npsva(
                self.args,
                variant,
                input_bam,
                sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )
            self.assertEqual(hom_alt_ref, 3.0)
            self.assertEqual(hom_alt_alt, 16.0)
            # Since this si a deletion, the total reference fragment count is divided by 2
            self.assertEqual(
                len(read_names["rl"]) + len(read_names["rr"]), hom_alt_ref * 2
            )
            self.assertEqual(len(read_names["al"]) + len(read_names["ar"]), hom_alt_alt)
            self.assertEqual(
                len(
                    set(read_names["rl"] + read_names["rr"])
                    & set(read_names["al"] + read_names["ar"])
                ),
                0,
            )

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_1.bam")
            sample = Sample.from_distribution(
                input_bam, 569.2, 163, 148, mean_coverage=25
            )
            het_ref, het_alt, *_ = count_alleles_with_npsva(
                self.args,
                variant,
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
            variant = Variant.from_pyvcf(record, None)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_distribution(
                input_bam, 569.2, 163, 148, mean_coverage=25
            )

            try:
                output_bam_file = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".bam"
                )
                output_bam_file.close()
                hom_alt_ref, hom_alt_alt, *_ = count_alleles_with_npsva(
                    self.args,
                    variant,
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

    def test_insert_distribution(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), input_bam)
            count_alleles_with_npsva(
                self.args,
                variant,
                input_bam,
                sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
            )

            count_alleles_with_npsva(
                self.args,
                variant,
                input_bam,
                sample,
                input_fasta=self.input_fasta,
                ref_contig="1_2073761_2073846_DEL",
                alt_contig="1_2073761_2073846_DEL_alt",
                insert_hist=False,
            )


class NPSVAFragmentTest(unittest.TestCase):
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

    def test_baseline_straddle_counting(self):
        for record in vcf.Reader(self.vcf_file):
            self.assertTrue(record.is_sv)
            variant = Variant.from_pyvcf(record, None)

            input_bam = os.path.join(FILE_DIR, "1_2073761_2073846_DEL_2.bam")
            sample = Sample.from_npsv(os.path.join(FILE_DIR, "stats.json"), input_bam)

            fragments = npsva.RealignedFragments(
                sample.mean_insert_size,
                sample.std_insert_size,
                sample.insert_size_density().as_dict(),
                input_bam,
            )
            num_reads = fragments.gather_reads(variant.region_string(flank=self.args.flank))
            self.assertEqual(fragments.size(), 254)

            print(variant.event_length)
            pair_results = fragments.count_baseline_straddlers(
                variant.region_string(), self.args.flank, -variant.event_length, 1.5, 10,
            )
            self.assertAlmostEqual(pair_results["alt_weighted_count"], 13.96, places=1)
            self.assertAlmostEqual(pair_results["insert_lower"], 0.0, places=2)
            self.assertAlmostEqual(pair_results["insert_upper"] / pair_results["insert_count"], 0.16, places=2)
