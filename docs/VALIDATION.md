# Validation contract

## Bit identity versus numerical acceptance

The strict oracle avoids CSV round-tripping and compares every float as
`uint64`, corrected counts/outliers as `int32`, and all ordered titles. The
checked-in run in `results/bit_identity.json` passes with zero mismatches.
A live integration test also creates a new oracle with R and compares it in a
temporary directory, so the gate is not limited to a stored result.

The exact guarantee is pinned to Windows x64, R 4.5.2, sctransform 0.4.3,
glmGamPoi 1.22.0, matrixStats 1.5.0 and the Rtools45 floating-point behavior.
Cross-platform builds retain the numerical acceptance contract below until
they have their own binary oracles.

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

Run the stored oracle and then regenerate it live:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' `
  benchmarks/bit_reference.R benchmark_data build/bit_oracle
python benchmarks/bit_identity.py benchmark_data build/bit_oracle `
  results/bit_identity.json --enforce
$env:SCT_TEST_RSCRIPT = 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe'
python -m pytest -q tests/test_r_bit_identity.py
```
