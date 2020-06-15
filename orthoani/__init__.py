import contextlib
import io
import typing
import os
import multiprocessing.pool
from typing import Dict, List, Iterator, Tuple

import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as BlastN
from Bio.Blast.Applications import NcbimakeblastdbCommandline as MakeBlastDB

from ._utils import ExitStack, BlockIterator, temppath

__all__ = ["orthoani", "orthoani_pairwise"]
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = (
    __import__("pkg_resources")
    .resource_string(__name__, "_version.txt")
    .decode("utf-8")
    .strip()
)


@contextlib.contextmanager
def _database(reference: os.PathLike) -> Iterator[None]:
    """A context managing a database for the given reference genome.

    The database is created when the context is entered, and deleted when it
    is exited.

    Arguments:
        reference (`str`): The path to a FASTA file containing a chopped
            genome to build a database for.
    """
    try:
        MakeBlastDB(dbtype="nucl", input_file=os.fspath(reference))()
        yield
    finally:
        os.remove("{}.nhr".format(os.fspath(reference)))
        os.remove("{}.nin".format(os.fspath(reference)))
        os.remove("{}.nsq".format(os.fspath(reference)))


def _chop(record: SeqRecord, dest: os.PathLike, blocksize: int) -> None:
    """Chop a single record into blocks and write it to ``dest``.

    Arguments:
        record (`~Bio.SeqRecord.SeqRecord`): The record to chop.
        dest (`~os.PathLike`): The path to the file where to write the blocks.
        blocksize (`int`): The size of block, in nucleotides.

    """
    with open(dest, mode="w") as d:
        if len(record) < blocksize:
            record = record[:0]
            record.id = "{}_{}".format(record.id, 1)
            record.description = ""
            Bio.SeqIO.write(record, d, "fasta")
        for i, block in enumerate(BlockIterator(record, blocksize)):
            block.id = "{}_{}".format(block.id, i)
            block.description = ""
            Bio.SeqIO.write(block, d, "fasta")


def _hits(
    query: os.PathLike,
    reference: os.PathLike,
    blocksize: int = 1020,
    threads: int = multiprocessing.cpu_count(),
) -> Dict[Tuple[str, str], float]:
    """Compute the hits from ``query`` to ``reference``.
    """
    cmd = BlastN(
        task="blastn",
        query=os.fspath(query),
        db=os.fspath(reference),
        evalue=1e-5,
        dust="no",
        xdrop_gap=150,
        penalty=-1,
        reward=1,
        num_alignments=1,
        num_threads=4,
        outfmt=5,
    )
    output = io.StringIO(cmd()[0])

    identities = {}
    for record in NCBIXML.parse(output):
        if record.alignments:
            hsps = record.alignments[0].hsps
            if all(hsp.align_length > 0.35 * blocksize for hsp in hsps):
                q = record.query
                r = record.alignments[0].hit_def
                identities[q, r] = (
                    sum(hsp.identities / hsp.align_length for hsp in hsps) / len(hsps)
                )

    return identities


def _orthoani(
    query:  os.PathLike,
    reference:  os.PathLike,
    blocksize: int
) -> float:
    """Compute the OrthoANI score for two chopped sequences.

    Arguments:
        query (`~os.PathLike`): The path to the chopped query.
        reference (`~os.PathLike`): The path to the chopped reference.
        blocksize (`int`): The size of the blocks used for computation.

    Caution:
        This function expects that a BLASTn database has been created for both
        of its inputs, next to the source FASTA file.
    """
    forward = _hits(query, reference, blocksize=blocksize)
    backward = _hits(reference, query, blocksize=blocksize)
    hits = { k:v for k,v in forward if (v,k) in backward }

    ani = 0.0
    for hit_q, hit_r in hits.items():
        ani += (backward[hit_r, hit_q] + forward[hit_q, hit_r]) / 2
    if hits:
        ani /= len(hits)
    return ani


def orthoani(
    reference: SeqRecord,
    query: SeqRecord,
    blocksize: int = 1020,
) -> float:
    with ExitStack() as ctx:
        # make the chopped file and the database for the first sequence
        chopped_r = ctx << temppath(suffix=".fa")
        _chop(reference, chopped_r, blocksize=blocksize)
        ctx << _database(chopped_r)
        # make the chopped file and the database for the first sequence
        chopped_q = ctx << temppath(suffix=".fa")
        _chop(query, chopped_q, blocksize=blocksize)
        ctx << _database(chopped_q)
        # return the orthoani score
        return _orthoani(chopped_q, chopped_r, blocksize)


def orthoani_pairwise(
    genomes: List[SeqRecord],
    blocksize: int = 1020
) -> Dict[Tuple[str, str], float]:
    with ExitStack() as ctx:
        # chop all inputs and make databases to avoid redoing it for each pair
        chopped = {}
        for genome in genomes:
            chopped[genome.id] = ctx << temppath(suffix=".fa")
            _chop(genome, chopped[genome.id], blocksize=blocksize)
            ctx << _database(chopped[genome.id])
        # query the results for each pair (using symmetricity)
        results = {}
        for i, g1 in enumerate(genomes):
            results[g1.id, g1.id] = 1.0
            for g2 in tqdm.tqdm(genomes[i+1:], leave=False):
                ani = _orthoani(chopped[g1.id], chopped[g2.id], blocksize=blocksize)
                results[g1.id, g2.id] = results[g2.id, g1.id] = ani
        # return the orthoani score for each pair
        return results
