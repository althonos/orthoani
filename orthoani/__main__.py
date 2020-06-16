"""Simple CLI interface to the OrthoANI algorithm.
"""
import argparse
import multiprocessing
import signal
import sys

import better_exceptions
from Bio.SeqIO import read

from . import orthoani, __name__, __version__

parser = argparse.ArgumentParser(
    prog="orthoani",
    description="Compute OrthoANI between two sequences in FASTA format"
)
parser.add_argument(
    '--traceback',
    action='store_true',
    help="display a complete traceback on program error",
)
parser.add_argument(
    '-V',
    '--version', 
    action='version', 
    version='{} {}'.format(__name__, __version__)
)
parser.add_argument(
    "-j", 
    "--jobs",
    help="the number of threads to use for BLASTn", 
    default=multiprocessing.cpu_count(),
    type=int,
)
parser.add_argument(
    "-q", 
    "--query",
    metavar="Q", 
    help="the path to the first sequence to process",
    required=True
)
parser.add_argument(
    "-r", 
    "--reference", 
    metavar="R",
    help="the path to the second sequence to process",
    required=True
)

args = parser.parse_args()
better_exceptions.hook()

try:
    query = read(args.query, "fasta")
    reference = read(args.reference, "fasta")
    print(orthoani(query, reference))
except KeyboardInterrupt:
    print("Interrupted")
    sys.exit(-signal.SIGINT)
except Exception as e:
    print(e)
    try:
        if args.traceback:
            raise
    finally:
        sys.exit(getattr(e, "errno", 1))
