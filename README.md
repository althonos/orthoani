# OrthoANI

*A Python implementation of the [OrthoANI](https://doi.org/10.1099/ijsem.0.000760) algorithm for nucleotide identity measurement.*

[![TravisCI](https://img.shields.io/travis/althonos/orthoani/master.svg?logo=travis&maxAge=600&style=flat-square)](https://travis-ci.org/althonos/orthoani/branches)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/orthoani/)
[![PyPI](https://img.shields.io/pypi/v/orthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/orthoani)
[![Wheel](https://img.shields.io/pypi/wheel/orthoani.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/orthoani/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/orthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/orthoani/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/orthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/orthoani/#files)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/orthoani/blob/master/CHANGELOG.md)
[![GitHub issues](https://img.shields.io/github/issues/althonos/orthoani.svg?style=flat-square&maxAge=600)](https://github.com/althonos/orthoani/issues)

## 🗺️ Overview

OrthoANI is a metric proposed by [Lee *et al.*](https://doi.org/10.1099/ijsem.0.000760)
in 2016 to improve computation of Average Nucleotide Identity. It uses
[BLASTn](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) to find orthologous
blocks in a pair of sequences, and then computes the average identity only
considering alignments of reciprocal orthologs.

![Algorithm](https://www.microbiologyresearch.org/docserver/fulltext/ijsem/66/2/000760-f1.gif)

This project is a reimplementation of the closed-source Java implementation
provided by the authors on [`ezbiocloud.net`](https://www.ezbiocloud.net/sw/oat).
It relies on [Biopython](https://biopython.org/) to handle the I/O and the


## 🔧 Installing

Installing with `pip` is the easiest:
```console
$ pip install orthoani
```

`orthoani` also requires the BLAST+ binaries to be installed on your machine
and available somewhere in your `$PATH`.

## 💡 Example

Use Biopython to load two FASTA files, and then `orthoani.orthoani` to compute
the OrthoANI metric between them:
```python
import orthoani
from Bio.SeqIO import read

genome_1 = read("sequence1.fa", "fasta")
genome_2 = read("sequence2.fa", "fasta")

ani = orthoani.orthoani(genome_1, genome_2)
```

## 🐏 Memory

`orthoani` uses the machine temporary folder to handle BLAST+ input and output
files, which is configurable through
[`tempfile.tempdir`](https://docs.python.org/3/library/tempfile.html#tempfile.tempdir).
On some systems (like ArchLinux), this filesystem can reside only in memory,
which means that even


## 📜 About

This library is provided under the open-source
[MIT license](https://choosealicense.com/licenses/mit/).
Please cite this library if you are using it in a scientific
context using the following DOI:
[**10.5281/zenodo.595572**](https://doi.org/10.5281/zenodo.595572)

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original OrthoANI authors. It was developed by
[Martin Larralde](https://github.com/althonos/orthoani) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/).*
