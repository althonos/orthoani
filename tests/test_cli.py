import os
import glob
import pathlib
import unittest
import subprocess
import sys
from subprocess import PIPE

from Bio.SeqIO import read

import orthoani


class TestCli(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pathlib.Path(__file__).parent / "data"
        cls.p1, cls.p2 = map(os.fspath, cls.data.glob("*.fna"))

    def test_orthoani_cli(self):
        # check the score we get is the same as the OrthoANI Java implementation
        args = [sys.executable, "-m", "orthoani", "-q", self.p1, "-r", self.p2]
        proc = subprocess.run(args, stdout=PIPE, cwd=self.data.parent.parent)
        proc.check_returncode()
        self.assertEqual(proc.stdout, b"0.5725\n")
