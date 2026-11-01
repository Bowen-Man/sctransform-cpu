# Third-party provenance

The strict Windows backend is distributed under **GPL-3.0-or-later** and
contains or derives from the following free-software components:

- R 4.5.2 `src/library/stats/src/portsrc.f` and small utilities derived from
  `port.c`, R Core Team, GPL-2.0-or-later.
- R 4.5 Mathlib positive-argument routines used by the compatibility
  callbacks, R Core Team, GPL-2.0-or-later.
- glmGamPoi 1.22.0 dispersion formulas and operation order, Bioconductor
  commit `908cde9d595259b4a317343ddfbcf20036ef0895`, GPL-3.0; portions are
  identified upstream as adapted from DESeq2 under LGPL-3.0-or-later.
- sctransform 0.4.3 scalar Rcpp expressions and pipeline operation order,
  Satija Lab, GPL-3.0-or-later.
- matrixStats 1.5.0 `rowVars` reduction and refinement order, Henrik Bengtsson,
  Artistic-2.0.
- The pinned MinGW exponential compatibility fragment is derived from
  mingw-w64 CRT code under ZPL-2.1, as documented in its source file.

The bundled `sctransform_cpu/native/sct_port_r45.dll` was built with Rtools45
GCC/GFortran 14.3.0. Its SHA-256 is
`22396be5857c7aa81b90eb9b7bf52843f522323f6c2a39ca262831ce5823d257`.

The original copyright notices, license identifiers, and source-level
attribution are retained in `sctransform_cpu/native/`. See `LICENSE` for the
distribution license.
