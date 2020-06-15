"""Simple CLI interface to the OrthoANI algorithm.
"""
import argparse
import multiprocessing
from Bio.SeqIO import read

from . import orthoani

parser = argparse.ArgumentParser()
parser.add_argument("-j", "--jobs", default=multiprocessing.cpu_count())
parser.add_argument("-q", "--query", required=True)
parser.add_argument("-r", "--reference", required=True)
args = parser.parse_args()

query = read(args.query, "fasta")
reference = read(args.reference, "fasta")
print(orthoani(query, reference))
