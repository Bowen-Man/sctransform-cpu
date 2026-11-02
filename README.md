# sctransform-cpu

An installable, R-free-runtime implementation of the validated SCTransform v2
CPU paths. This repository records the **preliminary reference work** for a
separate native Python CPU/GPU implementation: source auditing, exact oracle
construction, numerical edge cases, packaging and the first CPU performance
evidence.

It is also usable software today. The v0.4.0 wheels run the documented SCT v2
contract without an R installation. On the checked 800 x 240 PBMC3K-derived
fixture, Python compute was 6.88x faster than the pinned Windows R reference
and 15.46x faster in the pinned Linux run, with about 10--12x lower sampled
peak RSS. These are fixture- and machine-specific measurements, not universal
speedup claims; see the [benchmark report](docs/BENCHMARK_REPORT.md).

The public API and pipeline orchestration are Python/NumPy/SciPy. Bundled
Windows x64 (`.dll`) and Linux x86_64 (`.so`) compatibility libraries
reproduce the platform-pinned R/Rcpp, glmGamPoi, matrixStats and PORT
arithmetic where ordinary NumPy/SciPy evaluation differs by a few bits.
**R is not a runtime dependency.** The package is a native-backed Python
implementation, not an all-Python numerical reimplementation.

## Repository role

This GPL repository remains the exact CPU reference, live-R oracle harness and
reproducible baseline. It will receive correctness, compatibility, packaging
and security maintenance, but it is not the production home for the future
MIT-compatible AnnData/CuPy/RAPIDS backend. That implementation will be
developed independently from the published mathematics and black-box
contracts; no R/Fortran runtime or build dependency from this package will be
carried into it. See the [reference status](docs/REFERENCE_STATUS.md) and
[maintenance roadmap](docs/ROADMAP.md).

## Validated scope

Validated reference environments use R 4.5.2, sctransform 0.4.3,
glmGamPoi 1.22.0 and matrixStats 1.5.0. The two distinct native contracts are
Windows x64/Rtools45 and Ubuntu 26.04 x86_64/glibc 2.43/GCC 15.2.

- direct v2 `vst` path with the default fixed `log_umi` slope;
- additional numeric regularized covariates;
- one categorical batch variable with batch-specific coefficients;
- fixed-theta non-regularized covariates when every fitted theta is finite;
- initial and regularized model parameters;
- Pearson residuals and residual variance;
- R 4.5-compatible implicit cell/gene sampling and reusable title manifests;
- Poisson-boundary, rare-gene, clipping and extreme-library adversarial cases;
- Top-100/200/300 residual-variance feature sets;
- corrected UMI counts for default, covariate, batch and supported nonreg designs;
- a Python-native adapter for the targeted official `vst` return names;
- a Windows-pinned 128-case fixed-mean dispersion stress suite.

This is a version- and contract-pinned alpha release, not yet a claim of parity
for every SCT/Seurat workflow. In particular, official v2 itself cannot fit a
non-regularized stage when regularization produces `theta=Inf`; this package
fails closed for the same case. Alternative residual and fitting methods,
arbitrary R contrasts, macOS, other Linux distributions/libm versions and the
complete Seurat object lifecycle are not yet covered.

The pinned 800 x 240 regression fixture passes the strict binary gate on both
validated platforms: every compared binary64 value, integer output and
ordered title is identical to a fresh same-platform R oracle. This is a
deliberately narrow platform/version guarantee, not a claim that arbitrary
R/package/platform combinations are bit-identical. See [Linux validation](docs/LINUX_VALIDATION.md)
and [development rules](docs/DEVELOPMENT.md).

## Install

The v0.4.0 release wheels support the complete validated contract listed
above on Windows x64 and Linux x86_64. Install the wheel for your platform
directly from that versioned GitHub release:

```powershell
python -m pip install "https://github.com/Bowen-Man/sctransform-cpu/releases/download/v0.4.0/sctransform_cpu-0.4.0-py3-none-win_amd64.whl"
```

```bash
python -m pip install "https://github.com/Bowen-Man/sctransform-cpu/releases/download/v0.4.0/sctransform_cpu-0.4.0-py3-none-linux_x86_64.whl"
```

The project is not published on PyPI, so `pip install sctransform-cpu` without
the release URL is not currently supported. R is needed only to regenerate the
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
    config=SCTConfig(output="hvg", features=3000, random_seed=42),
)

# Genes returned as rows of result.residuals
residual_genes = result.genes[result.residual_gene_indices]
```

The v2 defaults use at most 2,000 fit cells and 2,000 fit genes, matching the
official flavor rewrite. When the input exceeds either limit, pass
`random_seed` to reproduce R 4.5's default sampling stream, or reuse an
explicit ordered selection manifest:

```python
from sctransform_cpu import SCTConfig, make_selection_manifest, sctransform

manifest = make_selection_manifest(
    counts,
    gene_names=gene_names,
    cell_names=cell_names,
    seed=42,
)
manifest.write("fit-selection.json")

result = sctransform(
    counts,
    gene_names=gene_names,
    cell_names=cell_names,
    config=SCTConfig(random_seed=42),
    selection_manifest=manifest,
)
```

The standalone `scripts/make_selection_manifest.py` command writes the same
ordered cell/gene titles without running GLM fitting. See the pinned
[v2 source map](docs/SCT_V2_SOURCE_MAP.md) for the exact filtering and RNG
contract.

Additional numeric covariates and one categorical batch variable use the same
`cell_attributes` mapping in cell order:

```python
# Numeric covariate: y ~ log_umi + percent_mito
numeric = sctransform(
    counts,
    cell_attributes={"percent_mito": percent_mito},
    config=SCTConfig(
        latent_variables=("log_umi", "percent_mito"),
        random_seed=42,
    ),
)

# Batch design: y ~ (log_umi):batch + batch + 0
batched = sctransform(
    counts,
    cell_attributes={"batch": batch_labels, "percent_mito": percent_mito},
    config=SCTConfig(
        batch_variable="batch",
        latent_variables_nonreg=("percent_mito",),
        return_corrected_umi=True,
        random_seed=42,
    ),
)
```

Covariates must be finite numeric vectors. Batch labels must be nonempty ASCII
strings, and a sampled fit-cell set must retain every level. General-design
coefficients are exposed by `coefficient_names`, `coefficients` and
`fit_coefficients`; the legacy `intercept`/`slope` convenience fields apply
only when those global columns exist. Non-regularized coefficients are exposed
separately by `nonreg_coefficient_names` and `nonreg_coefficients` because the
official algorithm does not smooth them. Corrected UMI uses the regularized
covariates at their medians while preserving observed batch membership.

The main result fields are:

| Field | Shape / type | Meaning |
|---|---|---|
| `genes`, `cells` | `(G,)`, `(C,)` strings | retained names in computation order |
| `theta` | `(G,)` float64 | regularized inverse overdispersion |
| `coefficient_names`, `coefficients` | names and `(G, P)` float64 | regularized design coefficients |
| `fit_coefficients` | `(F, P)` float64 | unregularized coefficients for fitted genes |
| `nonreg_coefficient_names`, `nonreg_coefficients` | names and optional `(G, Q)` float64 | fixed-theta post-regularization fit |
| `intercept`, `slope` | optional `(G,)` float64 | convenience views for global columns |
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
