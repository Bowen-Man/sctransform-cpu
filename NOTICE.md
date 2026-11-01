# Third-party provenance

The strict native backends are distributed under **GPL-3.0-or-later** and
contain or derive from the following free-software components:

- R 4.5.2 `src/library/stats/src/portsrc.f` and small utilities derived from
  `port.c`, R Core Team, GPL-2.0-or-later.
- R 4.5 Mathlib positive-argument routines used by the compatibility
  callbacks, R Core Team, GPL-2.0-or-later.
- glmGamPoi 1.22.0 dispersion formulas and operation order, Bioconductor
  commit `908cde9d595259b4a317343ddfbcf20036ef0895`, GPL-3.0; portions are
  identified upstream as adapted from DESeq2 under LGPL-3.0-or-later.
- glmGamPoi's ordered count-frequency table is reproduced with the platform
  C++ standard library because its `unordered_map` traversal order affects
  cancellation-sensitive dispersion estimates.
- sctransform 0.4.3 scalar Rcpp expressions and pipeline operation order,
  Satija Lab, GPL-3.0-or-later.
- matrixStats 1.5.0 `rowVars` reduction and refinement order, Henrik Bengtsson,
  Artistic-2.0.
- The pinned MinGW exponential compatibility fragment is derived from
  mingw-w64 CRT code under ZPL-2.1, as documented in its source file.

The bundled `sctransform_cpu/native/sct_port_r45.dll` was built with Rtools45
GCC/GFortran 14.3.0. Its SHA-256 is
`d62310f070a70891952eec610f95803a089e5a9cb38b2938a48cced38f74f04f`.

The bundled `sctransform_cpu/native/libsct_port_r45.so` was built on Ubuntu
26.04 x86_64 with GCC/GFortran 15.2.0 and glibc 2.43. Its SHA-256 is
`7ed38d800dbf70e97ff964b6b98795053cce5eef8f06990dd33f3e711a6570b1`.

The original copyright notices, license identifiers, and source-level
attribution are retained in `sctransform_cpu/native/`. See `LICENSE` for the
distribution license.
