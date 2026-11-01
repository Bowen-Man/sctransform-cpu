# Changelog

All notable changes are recorded here. The project follows semantic versioning
while its validated statistical contract remains explicitly version-pinned.

## 0.3.0 - 2026-11-02

- Ship the bit-pinned Linux x86_64 backend alongside Windows x64 support.
- Add strict same-platform R oracle evidence and Linux benchmark results.
- Add Windows/Linux functional CI for Python 3.10 through 3.13.
- Add isolated wheel smoke tests, native-source rebuild CI and tag-driven
  GitHub Release publishing with SHA-256 checksums.
- Expose package `__version__` and document the public result contract.

## 0.2.0 - 2026-10-31

- Match the pinned Windows R 4.5.2 SCT v2 reference bit-for-bit on the public
  regression fixture.

## 0.1.0 - 2026-10-29

- Publish the first validated CPU reference implementation.
