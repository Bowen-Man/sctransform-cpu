# sctransform-cpu

A compact, CPU-only Python implementation of the directly validated,
unbatched SCTransform v2 `vst` path.

The statistical pipeline is Python/NumPy/SciPy. A bundled Windows x64 library
reuses R 4.5.2's PORT state machine so flat Gamma–Poisson likelihoods follow
the same optimizer path as the pinned R reference. **R is not a runtime
dependency.**

## Validated scope

Reference environment: R 4.5.2, sctransform 0.4.3 and glmGamPoi 1.22.0.

- direct, one-batch v2 `vst` path;
- fixed `log_umi` slope and no additional covariates;
- initial and regularized model parameters;
- Pearson residuals and residual variance;
- Top-100/200/300 residual-variance feature sets;
- corrected UMI counts;
- 128 fixed-mean dispersion stress cases.

This is a version- and contract-pinned alpha release, not a claim of parity
for every SCT/Seurat workflow. Cross-language implicit sampling,
batch/covariates, other operating systems and the complete Seurat output
contract are not yet covered.

## Install and run

The first release supports Windows x64 and Python 3.12.

```powershell
python -m pip install .
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

On the checked-in 800 x 240 PBMC3K-derived fixture, the measured medians were:

| implementation | compute wall | process wall | sampled peak RSS |
|---|---:|---:|---:|
| Python | 0.858 s | 1.686 s | 75.5 MiB |
| R reference | 3.340 s | 4.451 s | 732.1 MiB |

All pointwise quality gates passed. These are small-fixture measurements, not
general speedup or memory claims. See the full [benchmark report](docs/BENCHMARK_REPORT.md).

## License and provenance

GPL-3.0-or-later. The source distribution contains the corresponding native
sources and build instructions. See [NOTICE.md](NOTICE.md).
