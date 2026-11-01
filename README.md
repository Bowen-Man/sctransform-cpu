# sctransform-cpu

A compact, CPU-only Python implementation of the directly validated,
unbatched SCTransform v2 `vst` path.

The public API and pipeline orchestration are Python/NumPy/SciPy. Bundled
Windows x64 (`.dll`) and Linux x86_64 (`.so`) compatibility libraries
reproduce the platform-pinned R/Rcpp, glmGamPoi, matrixStats and PORT
arithmetic where ordinary NumPy/SciPy evaluation differs by a few bits.
**R is not a runtime dependency.** The package is a native-backed Python
implementation, not an all-Python numerical reimplementation.

## Validated scope

Validated reference environments use R 4.5.2, sctransform 0.4.3,
glmGamPoi 1.22.0 and matrixStats 1.5.0. The two distinct native contracts are
Windows x64/Rtools45 and Ubuntu 26.04 x86_64/glibc 2.43/GCC 15.2.

- direct, one-batch v2 `vst` path;
- fixed `log_umi` slope and no additional covariates;
- initial and regularized model parameters;
- Pearson residuals and residual variance;
- Top-100/200/300 residual-variance feature sets;
- corrected UMI counts;
- a Windows-pinned 128-case fixed-mean dispersion stress suite.

This is a version- and contract-pinned alpha release, not a claim of parity
for every SCT/Seurat workflow. Cross-language implicit sampling,
batch/covariates, macOS, other Linux distributions/libm versions and the
complete Seurat output contract are not yet covered.

The pinned 800 x 240 regression fixture passes the strict binary gate on both
validated platforms: every compared binary64 value, integer output and
ordered title is identical to a fresh same-platform R oracle. This is a
deliberately narrow platform/version guarantee, not a claim that arbitrary
R/package/platform combinations are bit-identical. See [Linux validation](docs/LINUX_VALIDATION.md)
and [development rules](docs/DEVELOPMENT.md).

## Install and run

Validated wheels support Windows x64 and Linux x86_64. The strict release gate
currently uses Python 3.12; the runtime does not invoke R.

```powershell
python -m pip install .
```

```bash
python -m pip install ./sctransform_cpu-0.2.0-py3-none-linux_x86_64.whl
```

```python
from scipy import sparse
from sctransform_cpu import SCTConfig, sctransform

# raw_counts is genes x cells, with non-negative integer UMI counts
counts = sparse.csr_matrix(raw_counts)
result = sctransform(
    counts,
    gene_names=gene_names,
    cell_names=cell_names,
    config=SCTConfig(output="hvg", features=3000),
)
```

The default fits all eligible genes and cells. If subsampling is requested,
the package deliberately requires explicit ordered `fit_cell_indices` and
`fit_gene_names`; equal random seeds do not make R and Python sampling
identical on every input.

## Reproduce correctness and performance

```powershell
python -m pip install -e ".[test,benchmark]"
python -m pytest -q
python benchmarks/run_benchmark.py `
  --python python `
  --rscript "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" `
  --repeats 5
```

The benchmark performs one warm-up and five fresh-process measurements for
both Python and R. It records compute wall time, process wall time, CPU time
and sampled peak RSS. See [validation](docs/VALIDATION.md) and the checked-in
machine-readable reports under `results/`.

Linux uses the analogous command with `Rscript` and records its evidence under
`results/linux/`; R is needed only to regenerate the oracle or benchmark.

On the checked-in 800 x 240 PBMC3K-derived fixture, the measured medians were:

| implementation | compute wall | process wall | sampled peak RSS |
|---|---:|---:|---:|
| Python | 0.491 s | 1.195 s | 63.2 MiB |
| R reference | 3.380 s | 4.512 s | 731.3 MiB |

The Ubuntu 26.04 WSL2 run measured 0.206 s Python versus 3.185 s R compute
medians, with sampled peak RSS of 64.7 MiB versus 680.0 MiB. These numbers are
machine- and fixture-specific, not general speedup claims.

The strict binary and pointwise quality gates passed. These are small-fixture
measurements, not general speedup or memory claims. See the full
[benchmark report](docs/BENCHMARK_REPORT.md).

## License and provenance

GPL-3.0-or-later. The source distribution contains the corresponding native
sources and build instructions. See [NOTICE.md](NOTICE.md).
