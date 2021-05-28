"""Implementation of the OrthoANI algorithm for nucleotide identity measurement.
"""

import collections
import contextlib
import decimal
import io
import os
import shlex
import subprocess
import typing
from subprocess import PIPE, DEVNULL
from typing import Dict, List, Iterator, Iterable, Tuple, Union

try:
    from os import PathLike, fspath
except ImportError:  # python 3.5
    from pathlib import Path as PathLike  # type: ignore
    from builtins import str as fspath  # type: ignore

from ._utils import ExitStack, BlockIterator, temppath

if typing.TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord



__all__ = ["orthoani", "orthoani_pairwise"]
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__version__ = "0.5.0"


@contextlib.contextmanager
def _seqidlist(seqids: "Iterable[str]") -> "Iterator[PathLike[str]]":
    """Get a context manager to manage a sequence IDs file.

    The file is created when the context is entered, and deleted when it
    is exited.

    Arguments:
        seqids (iterable of `str`): The sequence IDs to write to the file.

    """
    with temppath(suffix=".pacc") as rawidfile_path:
        # write down sequence IDs
        with open(rawidfile_path, "w") as f:
            f.writelines(["{}\n".format(seqid) for seqid in seqids])

        # convert sequence IDs to the v5 BLASTDB format
        seqidfile_path = "{}.bsl".format(fspath(rawidfile_path))
        args = [
            "blastdb_aliastool",
            "-seqid_file_in",
            fspath(rawidfile_path),
            "-seqid_file_out",
            fspath(seqidfile_path),
        ]
        try:
            proc = subprocess.run(args, stdout=PIPE, stderr=PIPE)
            proc.check_returncode()
            yield seqidfile_path
        except subprocess.CalledProcessError as error:
            raise RuntimeError(proc.stderr or proc.stdout) from error
        finally:
            if os.path.exists(seqidfile_path):
                os.remove(seqidfile_path)


@contextlib.contextmanager
def _database(reference: "PathLike[str]") -> Iterator[None]:
    """Get a context to manage a database for the given reference genome.

    The database is created when the context is entered, and deleted when it
    is exited.

    Arguments:
        reference (`str`): The path to a FASTA file containing a chopped
            genome to build a database for.

    """
    args = [
        "makeblastdb",
        "-parse_seqids",
        "-dbtype",
        "nucl",
        "-in",
        fspath(reference),
        "-blastdb_version",
        "5",
    ]
    try:
        proc = subprocess.run(args, stderr=PIPE, stdout=DEVNULL)
        proc.check_returncode()
        yield
    except subprocess.CalledProcessError as error:
        raise RuntimeError(proc.stderr) from error
    finally:
        for ext in "nhr", "nin", "nsq", "ndb", "nog", "nos", "not", "ntf", "nto":
            path = os.path.extsep.join([fspath(reference), ext])
            if os.path.exists(path):
                os.remove(path)


def _chop(
    record: Union["SeqRecord", Iterable["SeqRecord"]],
    dest: "PathLike[str]",
    blocksize: int
) -> None:
    """Chop a single record into blocks and write it to ``dest``.

    Arguments:
        record (iterable of `~Bio.SeqRecord.SeqRecord`): The record to chop,
            wrapped in an iterable containing one or more contigs.
        dest (`~os.PathLike`): The path to the file where to write the blocks.
        blocksize (`int`): The size of block, in nucleotides.

    """
    if hasattr(record, "id") and hasattr(record, "seq"):
        record = [record]
    with open(fspath(dest), mode="w") as d:
        for r in record:
            for i, block in enumerate(BlockIterator(r, blocksize)):
                # skip blocks with more than 80% unknown nucleotides
                n_count = block.seq.count("N") + block.seq.count("n")
                if n_count / len(block) >= 0.8:
                    continue
                # write FASTA with 80 characters per line
                d.write(">{}_{}\n".format(r.id, i))
                for i in range(0, len(block), 80):
                    d.write(str(block.seq[i:i+80]))
                    d.write("\n")


def _hits(
    query: "PathLike[str]",
    reference: "PathLike[str]",
    blocksize: int,
    threads: int,
    seqids: "Optional[List[str]]" = None,
) -> Dict[Tuple[str, str], decimal.Decimal]:
    """Compute the hits from ``query`` to ``reference``."""
    # run BLASTn
    with ExitStack() as ctx:
        args = [
            'blastn',
            '-query', fspath(query),
            '-db', fspath(reference),
            '-task', 'blastn',
            '-evalue', '1e-15',
            '-xdrop_gap', '150',
            '-dust', 'no',
            '-penalty', '-1',
            '-reward', '1',
            '-num_alignments', '1',
            '-outfmt', '6 qseqid sseqid length pident',
            '-num_threads', str(threads),
        ]
        if seqids is not None:
            args.extend(["-seqidlist", ctx << _seqidlist(seqids)])

        try:
            proc = subprocess.run(args, stdout=PIPE, stderr=PIPE)
            proc.check_returncode()
        except subprocess.CalledProcessError as error:
            raise RuntimeError(proc.stderr or proc.stdout) from error

    # group alignments together by query/ref couple
    pidents = collections.defaultdict(list)
    for line in proc.stdout.splitlines():
        qseqid, sseqid, length, pident = line.decode().split("\t")
        if int(length) >= 0.35 * blocksize:
            pidents[qseqid, sseqid].append(decimal.Decimal(pident))

    # return all HSP identities
    return pidents


def _orthoani(
    query: "PathLike[str]", reference: "PathLike[str]", blocksize: int, threads: int
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
    # compute forward hits: every sequence is used
    forward = _hits(query, reference, blocksize=blocksize, threads=threads)

    # compute forward hits: only sequences that got a forward hit are used
    backward = _hits(
        reference,
        query,
        blocksize=blocksize,
        threads=threads,
        seqids={k for k, v in forward},
    )

    # find reciprocical hits
    hits = {(k, v) for k, v in forward if (v, k) in backward}

    ani, hsps = decimal.Decimal(0), 0
    for hit_q, hit_r in hits:
        ani += sum(backward[hit_r, hit_q]) + sum(forward[hit_q, hit_r])
        hsps += len(backward[hit_r, hit_q]) + len(forward[hit_q, hit_r])
    if hsps:
        ani /= hsps * 100
    return float(ani)


def orthoani(
    reference: Union["SeqRecord", Iterable["SeqRecord"]],
    query: Union["SeqRecord", Iterable["SeqRecord"]],
    blocksize: int = 1020,
    threads: int = os.cpu_count(),
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
        ani = _orthoani(chopped_q, chopped_r, blocksize, threads)
    return ani


def orthoani_pairwise(
    genomes: List["SeqRecord"],
    blocksize: int = 1020,
    threads: int = os.cpu_count(),
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
