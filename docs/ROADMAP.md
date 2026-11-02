# Reference repository status and maintenance roadmap

## Mission

This repository is the completed preliminary SCTransform v2 reference phase.
It provides an installable GPL CPU package, pinned Windows/Linux numerical
contracts, live-R oracle tooling and a reproducible small-workload benchmark.
It is useful independently and is also the correctness baseline for later
native Python CPU/GPU research.

The future MIT-compatible AnnData/CuPy/RAPIDS implementation belongs in a
separate repository. R, Fortran and GPL-derived native compatibility code from
this package must not become dependencies of that implementation.

## v0.4.0 completion state

| Area | State | Evidence or boundary |
|---|---|---|
| R-free installed runtime | complete | wheels do not import or invoke R |
| Windows/Linux packages | complete | platform wheels and isolated smoke tests |
| Default fixed-slope v2 path | validated | same-platform live R binary oracle |
| R-compatible cell/gene sampling | validated | ordered names and RNG decisions |
| Numeric covariates | validated | parameters, regularization and residuals |
| One categorical batch variable | validated | batch-specific fit and regularization |
| Finite-theta nonreg variables | validated | coefficients and final residuals |
| Corrected UMI | validated | default and supported general designs |
| Edge cases | validated in current contract | Poisson boundary, clipping, rare genes and extreme library |
| Small-fixture CPU evidence | complete | Windows/Linux timing, CPU and sampled RSS |
| Full historical SCT/Seurat API | out of scope | unsupported branches are rejected or omitted |
| AnnData/RAPIDS/GPU backend | transferred | separate MIT-compatible project |

## Maintenance policy

Accepted work in this repository:

- correctness fixes within the documented contract;
- new regression fixtures for a reported failure;
- compatibility builds with their own platform oracle;
- packaging, security, provenance and documentation maintenance;
- oracle or benchmark tooling needed by the separate implementation.

Not planned here:

- `nb_offset`, historical fit modes or every residual type;
- complete Seurat object lifecycle emulation;
- broad factor/contrast language compatibility;
- AnnData orchestration, CuPy, RMM or GPU kernels;
- large-scale execution redesign and async scheduling.

Those exclusions are explicit scope decisions, not claims that the features
exist. The detailed option-level map remains in `SCT_V2_SOURCE_MAP.md`.

## Evidence boundary

The exact guarantee is version and platform pinned. The public 800 x 240
fixture proves regression behavior for the covered paths; it does not prove
every possible count matrix, platform or package version. Performance ratios
are likewise measurements of the named machines and fixture, not general
speedup promises.

No GPU performance claim is made by this repository.
