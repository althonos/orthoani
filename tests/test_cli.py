import os
import glob
import pathlib
import unittest
import subprocess
import sys

from Bio.SeqIO import read

import orthoani


class TestCli(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pathlib.Path(__file__).parent / "data"
        cls.p1, cls.p2 = map(os.fspath, cls.data.glob("1852379.*.fna"))

    def shell(self, args):
        return subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            cwd=self.data.parent.parent,
        )

    def test_orthoani_cli(self):
        # check the score we get is the same as the OrthoANI Java implementation
        args = [sys.executable, "-m", "orthoani", "-q", self.p1, "-r", self.p2]
        proc = self.shell(args)
        proc.check_returncode()
        self.assertEqual(proc.stdout, b"0.5725\n")

    def test_return_code(self):
        args = [sys.executable, "-m", "orthoani", "-q", self.data, "-r", self.p2]
        proc = self.shell(args)
        self.assertEqual(proc.returncode, 21) # IsADirectoryError

        args = [sys.executable, "-m", "orthoani", "-q", "_", "-r", self.p2]
        proc = self.shell(args)
        self.assertEqual(proc.returncode, 2) # FileNotFoundError

    def test_return_code_traceback(self):
        args = [sys.executable, "-m", "orthoani", "-T", "-q", self.data, "-r", self.p2]
        proc = self.shell(args)
        self.assertEqual(proc.returncode, 21) # IsADirectoryError

        args = [sys.executable, "-m", "orthoani", "-T", "-q", "_", "-r", self.p2]
        proc = self.shell(args)
        self.assertEqual(proc.returncode, 2) # FileNotFoundError
