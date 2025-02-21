# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://github.com/althonos/pyorthoani/compare/v0.6.0...HEAD

## [v0.7.0] - 2025-02-21
[v0.7.0]: https://github.com/althonos/pyorthoani/compare/v0.5.0...v0.6.0
### Changed
- Rename `orthoani` to `pyorthoani`.
- Use `pyproject.toml` instead of `setup.py` to handle project with `setuptools`.

## [v0.6.0] - 2024-11-06
[v0.6.0]: https://github.com/althonos/pyorthoani/compare/v0.5.0...v0.6.0
### Changed
- Make `better-exceptions` an optional dependency.
### Fixed
- Discrepancy with OrthoANI values for high similarity genomes ([#2](https://github.com/althonos/orthoani/issues)).
- Remaining files of temporary BLAST databases not being removed.
### Removed
- Support for Python 3.5.

## [v0.5.0] - 2021-05-28
[v0.5.0]: https://github.com/althonos/pyorthoani/compare/v0.4.0...v0.5.0
### Fixed
- Exception messages not properly rendering with `--traceback` enabled.
- Use `os.cpu_count` instead of `multiprocessing.cpu_count` where applicable.
- Make `BlockIterator` zero-copy using indices instead of slices.
### Changed
- Use `seqidlist` to reduce number of blocks compared in backward search.
- Switch to use BLASTdb v5 (instead of v4 previously).
- Prefix temporary files with `orthoani` prefix.
- Skip blocks containing more than 80% unknown nucleotides like in original implementation.
- Average ANI values based on HSPs instead of reciprocical blocks like in original implementation.

## [v0.4.0] - 2020-06-26
[v0.4.0]: https://github.com/althonos/pyorthoani/compare/v0.3.2...v0.4.0
### Changed
- Entire sequences smaller than the given blocksize will be padded with *N*.
- `biopython` requirement was relaxed to `v1.73`.
### Fixed
- Code using builtin API not available in Python3.5.

## [v0.3.2] - 2020-06-24
[v0.3.2]: https://github.com/althonos/pyorthoani/compare/v0.3.1...v0.3.2
### Fixed
- `orthoani` CLI ignoring the value of the `--jobs` flag.

## [v0.3.1] - 2020-06-24
[v0.3.1]: https://github.com/althonos/pyorthoani/compare/v0.3.0...v0.3.1
### Fixed
- `blastn` and `makeblastdb` being called with `shell=True`, causing issues
  if shell cannot be forked.
- Traceback not being displayed even with `--traceback` flag.
### Changed
- ANI values are collected using `decimal.Decimal` instead of `float`.

## [v0.3.0] - 2020-06-22
[v0.3.0]: https://github.com/althonos/pyorthoani/compare/v0.2.1...v0.3.0
### Added
- Support for genomes segmented in multiple contigs.

## [v0.2.1] - 2020-06-19
[v0.2.1]: https://github.com/althonos/pyorthoani/compare/v0.2.0...v0.2.1
### Fixed
- Temporary files with chopped FASTA not being deleted.

## [v0.2.0] - 2020-06-16
[v0.2.0]: https://github.com/althonos/pyorthoani/compare/v0.1.0-post1...v0.2.0
### Added
- `threads` argument controlling BLASTn thread count to `orthoani` and `orthoani_pairwise`.
- `-j` / `--jobs` flag controlling BLASTn thread count to CLI.
- Proper documentation and error codes to CLI.
- `orthoani` console script to call the CLI without `python -m` invocation.

## [v0.1.0-post1] - 2020-06-15
[v0.1.0-post1]: https://github.com/althonos/pyorthoani/compare/v0.1.0...v0.1.0-post1
### Fixed
- Travis-CI badge not rendering on `README.md`.
### Changed
- Made development status *Beta* instead of *Alpha* in `setup.cfg`.

## [v0.1.0] - 2020-06-15
[v0.1.0]: https://github.com/althonos/pyorthoani/compare/21725fe...v0.1.9
Initial release.
