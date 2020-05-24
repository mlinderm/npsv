import unittest
from npsv import npsva

class BindingTest(unittest.TestCase):
    def test_add(self):
        # test that 1 + 1 = 2
        self.assertEqual(npsva.add(1, 1), 2)

if __name__ == '__main__':
    unittest.main()