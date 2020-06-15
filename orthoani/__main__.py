import argparse
from Bio.SeqIO import read

from . import orthoani


parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", required=True)
parser.add_argument("-r", "--reference", required=True)

args = parser.parse_args()

query = read(args.query, "fasta")
reference = read(args.reference, "fasta")
print(orthoani(query, reference))
