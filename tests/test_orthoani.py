import os
import glob
import pathlib
import unittest

from Bio.SeqIO import read

import orthoani


class TestOrthoani(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pathlib.Path(__file__).parent / "data"
        cls.r1, cls.r2 = map(lambda p: read(p, "fasta"), cls.data.glob("*.fna"))

    def test_orthoani_single(self):
        # check the score we get is the same as the OrthoANI Java implementation
        ani = orthoani.orthoani(self.r1, self.r2)
        self.assertAlmostEqual(ani, 0.5725)

    def test_orthoani_pairwise(self):
        mapping = orthoani.orthoani_pairwise([self.r1, self.r2])
        self.assertAlmostEqual(mapping[self.r1.id, self.r2.id], 0.5725)
        self.assertAlmostEqual(mapping[self.r2.id, self.r1.id], 0.5725)
        self.assertAlmostEqual(mapping[self.r1.id, self.r1.id], 1.0)
        self.assertAlmostEqual(mapping[self.r2.id, self.r2.id], 1.0)
