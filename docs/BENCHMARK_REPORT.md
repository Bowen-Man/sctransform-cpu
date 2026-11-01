# CPU comparison with the R reference

## Contract

Both implementations received the same 800-gene x 240-cell PBMC3K-derived
raw UMI matrix (`nnz=22,131`). Both retained 800 genes, fitted the same 441
eligible genes over all 240 cells, materialized the full residual matrix and
returned corrected UMI counts.

Runs used one warm-up followed by five measured fresh subprocesses. Python and
R order alternated between repetitions. `OMP_NUM_THREADS`, `MKL_NUM_THREADS`
and `OPENBLAS_NUM_THREADS` were all one. Compute time excludes input loading
and result serialization; process wall includes interpreter and package
startup. RSS was sampled every 10 ms and can miss shorter peaks.

## Performance

| implementation | compute wall median [min, max] | process wall median | CPU median | sampled peak RSS median |
|---|---:|---:|---:|---:|
| Python | 0.491 [0.486, 0.497] s | 1.195 s | 0.484 s | 63.2 MiB |
| R reference | 3.380 [3.370, 3.410] s | 4.512 s | 3.390 s | 731.3 MiB |

On this one small fixture, Python compute was 6.88x faster, fresh-process wall
was 3.78x faster, and sampled peak RSS was 11.58x lower. These ratios are not
claims about larger datasets or different machines. In particular, R package
startup and the richer R return object contribute to process wall and RSS.

## Numerical quality

All predefined pointwise gates passed:

| output | max absolute difference | tolerance failures |
|---|---:|---:|
| step-1 intercept | 0 | 0 |
| step-1 theta | 0 | 0 |
| regularized intercept | 0 | 0 |
| regularized theta | 0 | 0 |
| 192,000 Pearson residuals | 0 | 0 |
| residual variance | 0 | 0 |

The outlier mask had zero mismatches, corrected UMI counts were exactly equal,
the complete residual-variance order was identical and residual-variance
Spearman correlation was 1.0. The separate raw-binary oracle also reported
zero bit mismatches for every compared intermediate and output.

## Profile

Deterministic profiling perturbed Python compute from about 0.49 s to 0.914 s,
so profile times must not be mixed into the benchmark. Within the profiled
pipeline, `fit_batch` consumed 0.886 s cumulative and the 441 PORT dispersion
fits consumed 0.756 s cumulative. Their nested likelihood objectives consumed
0.442 s cumulative; nested rows cannot be added together. The evidence still
identifies parameter fitting—not residual generation (~0.008 s under
profiling)—as the first optimization target for this fixture.

Machine-readable evidence is in `results/benchmark.json`,
`results/quality.json`, `results/python_profile_top.json` and
`results/environment.json`.

## Linux replication

The same fixture and benchmark protocol were run in Ubuntu 26.04 WSL2 on an
AMD Ryzen 9 9950X3D, Python 3.12.15 and the pinned Linux R stack. All thread
environment variables remained one.

| implementation | compute wall median [min, max] | process wall median | CPU median | sampled peak RSS median |
|---|---:|---:|---:|---:|
| Python | 0.206 [0.195, 0.207] s | 0.365 s | 0.206 s | 64.7 MiB |
| R reference | 3.185 [3.107, 3.241] s | 4.216 s | 3.180 s | 680.0 MiB |

The Linux quality and raw-binary gates had zero differences. The stage timer
attributed 0.177 s to fitting and 0.0073 s to residual statistics in its one
quality run. These values reinforce the small-fixture hotspot observation but
do not establish scaling behavior. Evidence is under `results/linux/`; see
`docs/LINUX_VALIDATION.md` for the exact native and wheel contract.
