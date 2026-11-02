# Linux native contract

Linux uses a bundled ELF shared object, `libsct_port_r45.so`, in the same role
as the Windows compatibility DLL. It contains the R PORT state machine and the
narrow scalar/reduction routines needed to reproduce the official Linux
execution order. It does not embed, link to or invoke an R runtime.

## Pinned environment

- Ubuntu 26.04 LTS under WSL2, x86_64;
- Linux 6.18.33.2, glibc 2.43;
- AMD Ryzen 9 9950X3D, 16 physical / 32 logical cores;
- Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1;
- GCC/GFortran 15.2.0;
- R 4.5.2, sctransform 0.4.3, glmGamPoi 1.22.0, matrixStats 1.5.0;
- `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=1` for benchmarks.

The `.so` dynamically uses the distribution `libgfortran.so.5`, `libgcc_s`,
`libm` and `libc`; libstdc++ plus BLAS/LAPACK are statically linked. Its
SHA-256 is
`97b426f7beced380fd98ded949cd1ad736decff503638d29777616efaf1b3645`.

## Exact gate

The 800-gene x 240-cell fixture retained 800 genes and fitted 441 genes over
all 240 cells. A fresh Linux R oracle and Python agreed bit-for-bit for all
compared intermediates and outputs, including 192,000 values each for model
mean, model variance, numerator, denominator and Pearson residual; 192,000
corrected UMI integers; all parameters, residual statistics, masks, names and
the complete HVG order. The result is in
`results/linux/bit_identity.json`.

The same environment also passes the numeric-covariate and categorical-batch
live gates for fitted/regularized parameters, outlier masks, full Pearson
residuals and residual statistics.

The expanded gate additionally passes numeric-covariate corrected UMI and a
batch plus non-regularized covariate design. Nonreg coefficients, final
residuals, residual statistics and corrected UMI integers are bit-identical.

This is a regression guarantee for the named fixture and pinned environment,
not a proof for every input or another libc/libm. The independent Windows
oracle must not be reused on Linux because the official R packages themselves
take platform-specific scalar paths.

## Build, test and package

```bash
./scripts/build_native.sh
SCT_TEST_RSCRIPT=/usr/bin/Rscript python -m pytest -q
python -m pip install ".[release]"
python -m build --sdist --wheel --outdir release_dist
python scripts/smoke_wheel.py release_dist
```

The clean wheel is tagged `py3-none-linux_x86_64`, contains the `.so` but not
the Windows DLL, and was installed into a new Python 3.12 environment. The
installed wheel passed the same raw-binary oracle with zero mismatches and
found the shared object from `site-packages` without an override variable.

## Benchmark

One warm-up and five alternating fresh-process measurements on the small
fixture produced these medians:

| implementation | compute wall | process wall | sampled peak RSS |
|---|---:|---:|---:|
| Python wheel/source-equivalent path | 0.206 s | 0.365 s | 64.7 MiB |
| R reference | 3.185 s | 4.216 s | 680.0 MiB |

The single quality run spent 0.177 s of 0.192 s recorded stage time in model
fitting; residual statistics used 0.0073 s. The fixture is too small to infer
large-data scaling or a general speedup. Full runs and quality metrics are in
`results/linux/benchmark.json` and `results/linux/quality.json`.
