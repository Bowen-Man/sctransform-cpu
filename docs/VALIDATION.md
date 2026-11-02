# Validation contract

## Bit identity versus numerical acceptance

The strict oracle avoids CSV round-tripping and compares every float as
`uint64`, corrected counts/outliers as `int32`, and all ordered titles. The
checked-in run in `results/bit_identity.json` passes with zero mismatches.
A live integration test also creates a new oracle with R and compares it in a
temporary directory, so the gate is not limited to a stored result.

The exact guarantee has two separate platform contracts. Both pin R 4.5.2,
sctransform 0.4.3, glmGamPoi 1.22.0 and matrixStats 1.5.0. Windows pins the
Rtools45 floating-point behavior; Linux pins Ubuntu 26.04 x86_64, glibc 2.43
and GCC/GFortran 15.2. The Linux gate generates its own oracle and does not
compare against Windows bytes.

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

The Windows-pinned stored 128-case fixed-mean solver suite must have zero
tolerance failures and identical zero/nonzero alpha classification. It is not
reused as a Linux oracle. Each platform's end-to-end comparison requires
identical outlier masks and Top-100/200/300 feature sets.

The public benchmark fixture is an 800-gene by 240-cell raw UMI subset derived
from the public 10x PBMC3K dataset. It is a computational regression fixture,
not a biological validation dataset.

The exact contract covers the direct fixed-slope path, regularized numeric
covariates and one categorical batch variable. Numeric covariates use the
official offset model; the batch path uses the official batch-specific
intercepts and slopes. Unsupported cases remain explicit: non-regularized
latent variables, corrected UMI for general designs, arbitrary R contrast
coding, complete Seurat post-processing, other R/package versions, macOS and
other Linux runtimes require separate implementations and gates.

Run the stored oracle and then regenerate it live:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' `
  benchmarks/bit_reference.R benchmark_data build/bit_oracle
python benchmarks/bit_identity.py benchmark_data build/bit_oracle `
  results/bit_identity.json --enforce
$env:SCT_TEST_RSCRIPT = 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe'
python -m pytest -q tests/test_r_bit_identity.py
python -m pytest -q tests/test_r_design_identity.py
```

Linux equivalent:

```bash
./scripts/build_native.sh
export SCT_TEST_RSCRIPT=/usr/bin/Rscript
python -m pytest -q tests/test_r_bit_identity.py
python -m pytest -q tests/test_r_design_identity.py
Rscript benchmarks/bit_reference.R benchmark_data build/bit_oracle_linux
python benchmarks/bit_identity.py benchmark_data build/bit_oracle_linux \
  results/linux/bit_identity.json --enforce
```

The oracle exports the raw `genes_log_gmean` expression used inside `vst()`.
It intentionally does not apply `log10(10^genes_log_gmean)` through the public
`gene_attr$gmean` field, because that extra reporting round trip can introduce
one unused ULP and is not the value consumed by fitting or regularization.
