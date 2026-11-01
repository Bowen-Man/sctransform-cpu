# Validation contract

## Bit identity versus numerical acceptance

Version 0.1.0 passed the numerical gates below; it did **not** pass exact
IEEE-754 bit identity. A new binary oracle avoids CSV round-tripping and
compares every float as `uint64`, corrected counts/outliers as `int32`, and
all ordered titles. The first run is recorded in `results/bit_identity.json`
and currently fails. See `DEVELOPMENT.md` for the mismatch breakdown and
remediation order.

The acceptance tolerances are explicit engineering gates, not universal
biological error bounds:

| Quantity | absolute tolerance | relative tolerance |
|---|---:|---:|
| intercept | 1e-7 | 1e-7 |
| theta | 1e-8 | 1e-5 |
| Pearson residual | 1e-6 | 1e-5 |
| residual variance | 1e-7 | 1e-5 |
| fixed slope | 1e-13 | 0 |
| corrected UMI | 0 | 0 |

The stored 128-case fixed-mean solver suite must have zero tolerance failures
and identical zero/nonzero alpha classification. End-to-end comparison also
requires identical outlier masks and Top-100/200/300 feature sets.

The public benchmark fixture is an 800-gene by 240-cell raw UMI subset derived
from the public 10x PBMC3K dataset. It is a computational regression fixture,
not a biological validation dataset.

Limitations remain explicit: only the direct unbatched path is covered;
implicit cross-language sampling, additional covariates, complete Seurat
post-processing, other R/package versions and non-Windows builds require
separate gates.
