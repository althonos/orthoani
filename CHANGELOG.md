# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://github.com/althonos/orthoani/compare/v0.2.0...HEAD

## [v0.2.0] - 2020-06-16
[v0.2.0]: https://github.com/althonos/orthoani/compare/v0.1.0-post1...v0.2.0
### Added
- `threads` argument controlling BLASTn thread count to `orthoani` and `orthoani_pairwise`.
- `-j` / `--jobs` flag controlling BLASTn thread count to CLI.
- Proper documentation and error codes to CLI.
- `orthoani` console script to call the CLI without `python -m` invocation.

## [v0.1.0-post1] - 2020-06-15
[v0.1.0-post1]: https://github.com/althonos/orthoani/compare/v0.1.0...v0.1.0-post1
### Fixed
- Travis-CI badge not rendering on `README.md`.
### Changed
- Made development status *Beta* instead of *Alpha* in `setup.cfg`.

## [v0.1.0] - 2020-06-15
[v0.1.0]: https://github.com/althonos/orthoani/compare/21725fe...v0.1.9
Initial release.
