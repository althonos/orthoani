import os
import glob
import pathlib
import unittest
import tempfile

try:
    from os import fspath
except ImportError:
    from builtins import str as fspath

from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import read, parse

import orthoani


class TestOrthoani(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pathlib.Path(__file__).parent / "data"
        cls.r1, cls.r2 = map(lambda p: read(fspath(p), "fasta"), cls.data.glob("1852379.*.fna"))

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


class TestChop(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = pathlib.Path(__file__).parent / "data"

    def test_single_contig(self):
        record = read(fspath(self.data / "U00096-3.fasta"), "fasta")
        with (self.data / "U00096-3.fasta.chopped.fasta").open() as f:
            expected = list(parse(f, "fasta"))
        with tempfile.NamedTemporaryFile(mode="rt", suffix=".fna") as tmp:
            orthoani._chop(record, tmp.name, 1020)
            actual = list(parse(tmp, "fasta"))
        for actual_record, expected_record in zip(actual, expected):
            self.assertEqual(actual_record.seq, expected_record.seq)

    @unittest.skip("we pad smaller individual records")
    def test_multiple_contig(self):
        record = parse(fspath(self.data / "NZ_AAEN01000029.fna"), "fasta")
        with (self.data / "NZ_AAEN01000029.fna.chopped.fasta").open() as f:
            expected = list(parse(f, "fasta"))
        with tempfile.NamedTemporaryFile(mode="rt", suffix=".fna") as tmp:
            orthoani._chop(record, tmp.name, 1020)
            actual = list(parse(tmp, "fasta"))
        for actual_record, expected_record in zip(actual, expected):
            self.assertEqual(actual_record.seq, expected_record.seq)
