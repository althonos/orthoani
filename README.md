# PyOrthoANI [![Stars](https://img.shields.io/github/stars/althonos/orthoani.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/orthoani/stargazers)

*A Python implementation of the [OrthoANI](https://doi.org/10.1099/ijsem.0.000760) algorithm for nucleotide identity measurement.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyorthoani/test.yml?branch=master&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyorthoani/actions)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyorthoani/)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyorthoani?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyorthoani/)
[![PyPI](https://img.shields.io/pypi/v/pyorthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyorthoani)
[![Wheel](https://img.shields.io/pypi/wheel/pyorthoani.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyorthoani/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyorthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyorthoani/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyorthoani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyorthoani/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyorthoani/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyorthoani.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyorthoani/issues)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyorthoani/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyorthoani?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyorthoani)


## 🗺️ Overview

OrthoANI is a metric proposed by Lee *et al.*[\[1\]](#ref1)
in 2016 to improve computation of Average Nucleotide Identity. It uses
[BLASTn](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) to find orthologous
blocks in a pair of sequences, and then computes the average identity only
considering alignments of reciprocal orthologs.

![Algorithm](https://www.microbiologyresearch.org/docserver/fulltext/ijsem/66/2/000760-f1.gif)

PyOrthoANI is a reimplementation of the closed-source Java implementation
provided by the authors on [`ezbiocloud.net`](https://www.ezbiocloud.net/sw/oat).
It relies on [Biopython](https://biopython.org/) to handle the I/O, and calls
the BLAST+ binaries using the `subprocess` module of the Python standard 
library.


## 🔧 Installing

Installing with `pip` is the easiest:
```console
$ pip install pyorthoani
```

PyOrthoANI also requires the BLAST+ binaries to be installed on your machine
and available somewhere in your `$PATH`.


## 💡 Example

Use Biopython to load two FASTA files, and then `orthoani.orthoani` to compute
the OrthoANI metric between them:
```python
import pyorthoani
from Bio.SeqIO import read

genome_1 = read("sequence1.fa", "fasta")
genome_2 = read("sequence2.fa", "fasta")

ani = pyorthoani.orthoani(genome_1, genome_2)
```

`pyorthoani` can also be used from the CLI using a very simple command-line
interface mimicking the original Java tool:
```console
$ pyorthoani -q sequence1.fa -r sequence2.fa
57.25
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


## 🔖 Citation

PyOrthoANI is scientific software; it is submitted for publication
and is currently available as a [pre-print on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.02.13.638148v1).
Please cite both [PyOrthoANI](https://www.biorxiv.org/content/10.1101/2025.02.13.638148v1)
and [OrthoANI](https://pubmed.ncbi.nlm.nih.gov/26585518/) if you are using it in an academic work,
for instance as:

> PyOrthoANI (Larralde *et al.*, 2024), a Python implementation of OrthoANI (Lee *et al.*, 2016).


## 📜 About

This library is provided under the open-source
[MIT license](https://choosealicense.com/licenses/mit/).

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the [original OrthoANI authors](http://www.chunlab.com/). It was developed by
[Martin Larralde](https://github.com/althonos/orthoani) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*

## 📚 References

- <a id="ref1">\[1\]</a> Imchang Lee, Yeong Ouk Kim, Sang-Cheol Park and Jongsik Chun. *OrthoANI: An improved algorithm and software for calculating average nucleotide identity* (2016). International Journal of Systematic and Evolutionary Microbiology. [doi:10.1099/ijsem.0.000760](https://doi.org/10.1099/ijsem.0.000760). [PMID:26585518](https://pubmed.ncbi.nlm.nih.gov/26585518/).

