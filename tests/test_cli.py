import io
import os
import glob
import pathlib
import unittest
import subprocess
import sys
from unittest import mock
from os import fspath

from Bio.SeqIO import read

import orthoani._main


class TestCli(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = pathlib.Path(__file__).parent / "data"
        cls.p1, cls.p2 = map(fspath, cls.data.glob("1852379.*.fna"))

    def setUp(self):
        self.stdout = io.StringIO()
        self._stdout_mock = mock.patch("sys.stdout", new=self.stdout)
        self._stdout_mock.__enter__()
        self.stderr= io.StringIO()
        self._stderr_mock = mock.patch("sys.stderr", new=self.stderr)
        self._stderr_mock.__enter__()

    def tearDown(self):
        self._stdout_mock.__exit__(None, None, None)
        self._stderr_mock.__exit__(None, None, None)

    def shell(self, args):
        return orthoani._main.main(args)

    def test_orthoani_cli(self):
        # check the score we get is the same as the OrthoANI Java implementation
        args = ["-q", self.p1, "-r", self.p2]
        retcode = self.shell(args)
        self.assertEqual(retcode, 0)
        self.assertEqual(self.stdout.getvalue(), "57.25\n")

    def test_return_code(self):
        args = ["-q", fspath(self.data), "-r", self.p2]
        retcode = self.shell(args)
        self.assertEqual(retcode, 21) # IsADirectoryError

        args = ["-q", "_", "-r", self.p2]
        retcode = self.shell(args)
        self.assertEqual(retcode, 2) # FileNotFoundError

    def test_return_code_traceback(self):
        args = ["-T", "-q", fspath(self.data), "-r", self.p2]
        retcode = self.shell(args)
        self.assertEqual(retcode, 21) # IsADirectoryError

        args = ["-T", "-q", "_", "-r", self.p2]
        retcode = self.shell(args)
        self.assertEqual(retcode, 2) # FileNotFoundError
