# Changelog

All notable changes are recorded here. The project follows semantic versioning
while its validated statistical contract remains explicitly version-pinned.

## Unreleased

## 0.4.0 - 2026-11-02

- Finalize this repository as the installable GPL CPU reference and
  preliminary performance study for a separate native CPU/GPU project.
- Add explicit phase gates and a compatibility roadmap.
- Add an R 4.5-compatible Mersenne-Twister/rejection sampler, inverse-density
  gene sampling, and reusable ordered selection manifests.
- Add the pinned source/option map and Python-native `vst` return adapter.
- Add live-R implicit-sampling and adversarial Poisson/clipping/extreme-count
  gates on Windows and Linux.
- Preserve glmGamPoi's native unordered count-table order, closing a
  collision-dependent dispersion mismatch found by the adversarial suite.
- Add regularized numeric covariates and one categorical batch variable with
  strict same-platform R gates for fitted parameters, regularization,
  outliers, Pearson residuals and residual statistics.
- Add finite-theta non-regularized GLMs and corrected UMI for supported
  general designs, with strict Windows/Linux R gates and explicit failure for
  the official `theta=Inf` incompatibility.

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
