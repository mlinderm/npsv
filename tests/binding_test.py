import os, unittest
from npsv import npsva

FILE_DIR = os.path.join(os.path.dirname(__file__), "data")

class NPSVABindingTest(unittest.TestCase):
    def test_count_alignments(self):
        ar = npsva.AlleleReference(
            os.path.join(FILE_DIR, "1_899922_899992_DEL_1.synth.fasta")
        )

        counts = ar.count_alignments(
            os.path.join(FILE_DIR, "1_899922_899992_DEL_1.synth.bam"),
            "1_899922_899992_DEL:3000-3001",
            "1_899922_899992_DEL_alt:3000-3001",
            rr_region="1_899922_899992_DEL:3070-3071",
        )
        self.assertEqual(
            counts, {"rl_reads": 17, "rr_reads": 29, "al_reads": 14, "ar_reads": 0}
        )

