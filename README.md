# OrthoANI [![Stars](https://img.shields.io/github/stars/althonos/orthoani.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/orthoani/stargazers)

*A Python implementation of the [OrthoANI](https://doi.org/10.1099/ijsem.0.000760) algorithm for nucleotide identity measurement.*

[![TravisCI](https://img.shields.io/travis/com/althonos/orthoani/master?logo=travis&maxAge=600&style=flat-square)](https://travis-ci.com/althonos/orthoani/branches)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/orthoani/)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/orthoani?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/orthoani/)
[![Sanity](https://img.shields.io/codacy/grade/4a427dadd1864c93ab9a55cb34c389a0.svg?style=flat-square&maxAge=3600)](https://codacy.com/app/althonos/orthoani)
[![PyPI](https://img.shields.io/pypi/v/orthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/orthoani)
[![Wheel](https://img.shields.io/pypi/wheel/orthoani.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/orthoani/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/orthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/orthoani/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/orthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/orthoani/#files)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/orthoani/blob/master/CHANGELOG.md)
[![GitHub issues](https://img.shields.io/github/issues/althonos/orthoani.svg?style=flat-square&maxAge=600)](https://github.com/althonos/orthoani/issues)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Forthoani)](https://pepy.tech/project/orthoani)


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

`orthoani` can also be used from the CLI using a very simple command-line
interface:
```console
$ orthoani -q sequence1.fa -r sequence2.fa
0.5725
```


## 🐏 Memory

`orthoani` uses the machine temporary folder to handle BLAST+ input and output
files, which is configurable through
[`tempfile.tempdir`](https://docs.python.org/3/library/tempfile.html#tempfile.tempdir).
On some systems (like ArchLinux), this filesystem can reside in memory, which means
that your computer could have trouble processing very large files. If this
happens, try changing the value of the `tempfile.tempdir` to a directory that
is actually located on physical storage.


## 📏 Precision

Values computed by this package and the original Java implementation may differ
slightly because in Java the authors perform rounding of floating-point values
at the sub-percent level, while this library uses the full values.


## 📜 About

This library is provided under the open-source
[MIT license](https://choosealicense.com/licenses/mit/).

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the [original OrthoANI authors](http://www.chunlab.com/). It was developed by
[Martin Larralde](https://github.com/althonos/orthoani) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
