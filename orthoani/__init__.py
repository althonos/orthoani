"""Implementation of the OrthoANI algorithm for nucleotide identity measurement.
"""

import collections
import contextlib
import decimal
import io
import os
import multiprocessing.pool
import shlex
import subprocess
from typing import Dict, List, Iterator, Iterable, Tuple, Union

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


def _is_atgc(char: str):
    n = ord(char)
    return n == 65 or n == 84 or n == 67 or n == 71


@contextlib.contextmanager
def _database(reference: os.PathLike) -> Iterator[None]:
    """Get a context to manage a database for the given reference genome.

    The database is created when the context is entered, and deleted when it
    is exited.

    Arguments:
        reference (`str`): The path to a FASTA file containing a chopped
            genome to build a database for.

    """
    try:
        cmd = MakeBlastDB(dbtype="nucl", input_file=os.fspath(reference))
        args = shlex.split(str(cmd))
        subprocess.run(args, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        yield
    finally:
        os.remove("{}.nhr".format(os.fspath(reference)))
        os.remove("{}.nin".format(os.fspath(reference)))
        os.remove("{}.nsq".format(os.fspath(reference)))


def _chop(record: Union[SeqRecord, Iterable[SeqRecord]], dest: os.PathLike, blocksize: int) -> None:
    """Chop a single record into blocks and write it to ``dest``.

    Arguments:
        record (iterable of `~Bio.SeqRecord.SeqRecord`): The record to chop,
            wrapped in an iterable containing one or more contigs.
        dest (`~os.PathLike`): The path to the file where to write the blocks.
        blocksize (`int`): The size of block, in nucleotides.

    """
    if isinstance(record, SeqRecord):
        record = [record]
    with open(dest, mode="w") as d:
        for r in record:
            for i, block in enumerate(BlockIterator(r, blocksize)):
                block.id = "{}_{}".format(block.id, i)
                block.description = ""
                Bio.SeqIO.write(block, d, "fasta")


def _hits(
    query: os.PathLike, reference: os.PathLike, blocksize: int, threads: int,
) -> Dict[Tuple[str, str], decimal.Decimal]:
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
        max_target_seqs=1,
        num_threads=threads,
        outfmt=5,
    )
    args = shlex.split(str(cmd))
    proc = subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    hits = {}
    for record in NCBIXML.parse(io.BytesIO(proc.stdout)):
        if record.alignments:
            nid = length = decimal.Decimal(0)
            for hsp in record.alignments[0].hsps:
                if hsp.align_length >= 0.35 * blocksize:
                    nid += hsp.identities
                    length += hsp.align_length
            if length > 0:
                hits[record.query, record.alignments[0].hit_def] = nid / length
    return hits


def _orthoani(
    query: os.PathLike, reference: os.PathLike, blocksize: int, threads: int
) -> float:
    """Compute the OrthoANI score for two chopped sequences.

    Arguments:
        query (`~os.PathLike`): The path to the chopped query.
        reference (`~os.PathLike`): The path to the chopped reference.
        blocksize (`int`): The size of the blocks used for computation.
        threads (`int`): The number of threads to use to run ``blastn``.

    Caution:
        This function expects that a BLASTn database has been created for
        both of its inputs, next to the source FASTA file.

    """
    forward = _hits(query, reference, blocksize=blocksize, threads=threads)
    backward = _hits(reference, query, blocksize=blocksize, threads=threads)
    hits = {k: v for k, v in forward if (v, k) in backward}

    ani = decimal.Decimal(0)
    for hit_q, hit_r in hits.items():
        ani += backward[hit_r, hit_q] + forward[hit_q, hit_r]
    if hits:
        ani /= len(hits) * 2
    return float(ani)


def orthoani(
    reference: Union[SeqRecord, Iterable[SeqRecord]],
    query: Union[SeqRecord, Iterable[SeqRecord]],
    blocksize: int = 1020,
    threads: int = multiprocessing.cpu_count(),
) -> float:
    """Compute the OrthoANI score for two sequence records.

    Arguments:
        reference (`~Bio.SeqRecord.SeqRecord`): The first record to process.
            If an iterable is given, it is assumed to yield contigs of the same
            genome.
        query (`~Bio.SeqRecord.SeqRecord`): The second record to process.
            If an iterable is given, it is assumed to yield contigs of the same
            genome.
        blocksize (`int`): The size of blocks to use for computation. Defaults
            to *1020bp*, the size used in the OrthoANI reference paper.
        threads (`int`): The number of threads to use to run ``blastn``.
            Defaults to number of CPUs.

    Returns:
        `float`: The OrthoANI score as a floating-point number between 0 and 1.

    Note:
        ``reference`` and ``query`` can be exchanged, since OrthoANI is a
        symmetric measurement. Names are only retained for consistency with
        BLAST but bear no semantic value.

    """
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
        return _orthoani(chopped_q, chopped_r, blocksize, threads)


def orthoani_pairwise(
    genomes: List[SeqRecord],
    blocksize: int = 1020,
    threads: int = multiprocessing.cpu_count(),
) -> Dict[Tuple[str, str], float]:
    """Compute pairwise OrthoANI scores for a list of sequence records.

    Use this function instead of `orthoani.orthoani` to avoid chopping the
    inputs and building the same BLAST database several times.

    Arguments:
        genomes (`list` of `~Bio.SeqRecord.SeqRecord`): The records to process.
        blocksize (`int`): The size of blocks to use for computation. Defaults
            to *1020bp*, the size used in the OrthoANI reference paper.
        threads (`int`): The number of threads to use to run ``blastn``.
            Defaults to number of CPUs.

    Returns:
        `dict`: A dictionary that maps the OrthoANI score to a pair of
        record identifiers.

    """
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
            for g2 in genomes[i + 1 :]:
                ani = _orthoani(chopped[g1.id], chopped[g2.id], blocksize, threads)
                results[g1.id, g2.id] = results[g2.id, g1.id] = ani
        # return the orthoani score for each pair
        return results
