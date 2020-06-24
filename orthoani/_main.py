"""Simple CLI interface to the OrthoANI algorithm.
"""
import argparse
import multiprocessing
import signal
import sys

import better_exceptions
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse

from . import orthoani, __name__, __version__


def argument_parser():
    parser = argparse.ArgumentParser(
        prog="orthoani",
        description="Compute OrthoANI between two sequences in FASTA format.",
        epilog="Full documentation <https://github.com/althonos/orthoani>",
    )
    parser.add_argument(
        "-T",
        "--traceback",
        action="store_true",
        help="display a complete traceback on program error",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version="{} {}".format(__name__, __version__),
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
        required=True,
    )
    parser.add_argument(
        "-r",
        "--reference",
        metavar="R",
        help="the path to the second sequence to process",
        required=True,
    )
    return parser


def main(argv=None):
    parser = argument_parser()
    args = parser.parse_args(argv)

    try:
        query = parse(args.query, "fasta")
        reference = parse(args.reference, "fasta")
        print(orthoani(query, reference))
        return 0

    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return -signal.SIGINT

    except Exception as e:
        if args.traceback: 
            print(
                better_exceptions.format_exception(type(e), e, e.__traceback__), 
                file=sys.stderr
            )
        else:
            print(e, file=sys.stderr)
        return getattr(e, "errno", 1)
