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
        library_object = sample_object.get_library("HG002")
        self.assertIsNotNone(library_object)

        self.assertEqual(library_object.read_length, 148)
        self.assertAlmostEqual(library_object.mean_insert_size, 573.060562, places=1)
        self.assertAlmostEqual(library_object.std_insert_size, 164.215239, places=1)

        self.assertEqual(
            library_object.insert_size_density[10000],
            norm.pdf(10000, loc=573.060562, scale=164.215239),
        )

        # Compute search distance
        self.assertGreater(sample_object.search_distance(percentile=0.99), norm.ppf(0.99, library_object.mean_insert_size, library_object.std_insert_size))

        # Compute the normalized GC coverage
        self.assertAlmostEqual(
            library_object.gc_normalized_coverage(0.40), 1.110604, places=3
        )
        self.assertAlmostEqual(
            library_object.gc_normalized_coverage(0.4012345), 1.110604, places=3
        )

        # Find the maximum GC normalization factor
        max_gc = sample_object.max_gc_normalized_coverage(limit=2.0)
        self.assertGreaterEqual(max_gc, 1.0)
        self.assertLessEqual(max_gc, 2.0)

        max_gc = sample_object.max_gc_normalized_coverage(limit=2.0, start=0.2, stop=0.8)
        print(max_gc)
