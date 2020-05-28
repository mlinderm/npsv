import os, unittest
from scipy.stats import norm
from npsv.sample import Sample

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")

class SampleLoadTestSuite(unittest.TestCase):
    def test_loads_npsv_json(self):
        json_path = os.path.join(FILE_DIR, "stats.json")
        bam_path = os.path.join(FILE_DIR, "1_1598414_1598580_DEL.bam")

        # Override the BAM path to be specific to test file
        sample_object = Sample.from_npsv(json_path, bam_path=bam_path)
        self.assertTrue(sample_object.has_read_group("NA12878"))
        self.assertEqual(sample_object.name, "HG002")
        self.assertAlmostEqual(sample_object.mean_coverage, 25.46, places=1)

        # Get generic library
        library_object = sample_object.get_lib("HG002")
        self.assertIsNotNone(library_object)

        self.assertEqual(library_object.read_length, 148)
        self.assertAlmostEqual(library_object.mean_insert_size, 572.32, places=1)
        self.assertAlmostEqual(library_object.std_insert_size, 156.41, places=1)

        self.assertEqual(
            library_object.insert_size_density[10000],
            norm.pdf(10000, loc=572.32, scale=156.41),
        )

        # Compute the normalized GC coverage
        self.assertAlmostEqual(
            library_object.gc_normalized_coverage(0.40), 1.110604, places=3
        )
        self.assertAlmostEqual(
            library_object.gc_normalized_coverage(0.4012345), 1.110604, places=3
        )

