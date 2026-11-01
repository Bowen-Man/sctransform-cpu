# sctransform-cpu

An R-free, Python-first implementation of SCTransform v2, built around a
strictly validated CPU reference. The long-term goal is a complete Python SCT
v2 implementation with an optional GPU backend that interoperates naturally
with AnnData, CuPy and the RAPIDS single-cell ecosystem.

The current release is the correctness-first CPU foundation. It does **not**
yet contain a GPU backend, and it covers a narrower contract than the complete
official R workflow. Correctness coverage and CPU performance must be closed
and measured before a GPU execution design is selected.

The public API and pipeline orchestration are Python/NumPy/SciPy. Bundled
Windows x64 (`.dll`) and Linux x86_64 (`.so`) compatibility libraries
reproduce the platform-pinned R/Rcpp, glmGamPoi, matrixStats and PORT
arithmetic where ordinary NumPy/SciPy evaluation differs by a few bits.
**R is not a runtime dependency.** The package is a native-backed Python
implementation, not an all-Python numerical reimplementation.

## Project direction

Development follows three gated phases:

1. **R alignment:** reproduce the targeted `sctransform::vst` v2 statistical
   paths, sampling decisions, outputs and edge cases against pinned live R
   oracles.
2. **CPU implementation and optimization:** establish scaling profiles, then
   optimize only measured hotspots while preserving the reference gates.
3. **RAPIDS-compatible GPU backend:** design the backend and dataflow from the
   completed CPU evidence, including transfer, sparse/dense and VRAM costs.

The GPU phase will not assume that CuPy, Triton, CUDA C++ or whole-pipeline GPU
execution is automatically best. See the explicit [project roadmap](docs/ROADMAP.md)
and current [development rules](docs/DEVELOPMENT.md).

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

This is a version- and contract-pinned alpha release, not yet a claim of parity
for every SCT/Seurat workflow. Cross-language implicit sampling,
batch/covariates, macOS, other Linux distributions/libm versions and the
complete Seurat output contract are not yet covered.

The pinned 800 x 240 regression fixture passes the strict binary gate on both
validated platforms: every compared binary64 value, integer output and
ordered title is identical to a fresh same-platform R oracle. This is a
deliberately narrow platform/version guarantee, not a claim that arbitrary
R/package/platform combinations are bit-identical. See [Linux validation](docs/LINUX_VALIDATION.md)
and [development rules](docs/DEVELOPMENT.md).

## Install

Release wheels support Windows x64 and Linux x86_64. Install the wheel for
your platform directly from the versioned GitHub release:

```powershell
python -m pip install "https://github.com/Bowen-Man/sctransform-cpu/releases/download/v0.3.0/sctransform_cpu-0.3.0-py3-none-win_amd64.whl"
```

```bash
python -m pip install "https://github.com/Bowen-Man/sctransform-cpu/releases/download/v0.3.0/sctransform_cpu-0.3.0-py3-none-linux_x86_64.whl"
```

For development, clone the repository and run `python -m pip install -e .`.
The project is not yet published on PyPI, so `pip install sctransform-cpu` is
not currently a supported command. R is needed only to regenerate the
correctness oracle; it is never imported or invoked at runtime.

| Platform | functional CI | strict bit-identity environment |
|---|---|---|
| Windows x64 | Python 3.10-3.13 | Python 3.12, Windows/Rtools45 contract |
| Linux x86_64 | Python 3.10-3.13 on Ubuntu 24.04 | Python 3.12, Ubuntu 26.04 contract |

Functional CI verifies installation and tests. The narrower Python 3.12 rows
are the environments in which a fresh R oracle has been compared bit-for-bit.

## Quick start

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

# Genes returned as rows of result.residuals
residual_genes = result.genes[result.residual_gene_indices]
```

The default fits all eligible genes and cells. If subsampling is requested,
the package deliberately requires explicit ordered `fit_cell_indices` and
`fit_gene_names`; equal random seeds do not make R and Python sampling
identical on every input.

The main result fields are:

| Field | Shape / type | Meaning |
|---|---|---|
| `genes`, `cells` | `(G,)`, `(C,)` strings | retained names in computation order |
| `intercept`, `theta`, `slope` | `(G,)` float64 | regularized v2 model parameters |
| `residual_mean`, `residual_variance` | `(G,)` float64 | per-gene Pearson residual statistics |
| `hvg_indices` | `(min(features, G),)` integers | stable residual-variance ranking |
| `residuals` | `(features, C)` or `(G, C)` float64 | dense clipped Pearson residuals |
| `residual_gene_indices` | one row index per residual row | maps residual rows back to `genes` |
| `corrected_umi` | optional CSR matrix | populated when `return_corrected_umi=True` |
| `timing`, `diagnostics` | dictionaries | stage timings and fit-outlier metadata |

Input must be a raw, non-negative integer UMI count matrix in **genes x
cells** orientation. Do not pass log-normalized expression. `output="hvg"`
avoids retaining all-gene residuals; `output="full"` returns the full dense
residual matrix and can require `8 * G * C` bytes before other temporaries.

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

Contributors should follow the [development rules](docs/DEVELOPMENT.md) and
[release checklist](docs/RELEASING.md). Release artifacts include source and
SHA-256 checksums.
