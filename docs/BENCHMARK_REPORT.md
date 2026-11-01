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
| Python | 0.858 [0.853, 0.860] s | 1.686 s | 0.859 s | 75.5 MiB |
| R reference | 3.340 [3.320, 3.380] s | 4.451 s | 3.320 s | 732.1 MiB |

On this one small fixture, Python compute was 3.89x faster, fresh-process wall
was 2.64x faster, and sampled peak RSS was 9.70x lower. These ratios are not
claims about larger datasets or different machines. In particular, R package
startup and the richer R return object contribute to process wall and RSS.

## Numerical quality

All predefined pointwise gates passed:

| output | max absolute difference | tolerance failures |
|---|---:|---:|
| step-1 intercept | 8.88e-16 | 0 |
| step-1 theta | 8.63e-5 | 0 |
| regularized intercept | 1.78e-14 | 0 |
| regularized theta | 1.48e-8 | 0 |
| 192,000 Pearson residuals | 3.46e-8 | 0 |
| residual variance | 1.96e-8 | 0 |

The outlier mask had zero mismatches, corrected UMI counts were exactly equal,
and Top-100/200/300 residual-variance feature sets were identical. Residual
variance Spearman correlation was 0.9999999941.

## Profile

Deterministic profiling perturbed Python compute from about 0.86 s to 2.03 s,
so profile times must not be mixed into the benchmark. Within the profiled
pipeline, `fit_batch` consumed 1.985 s cumulative and the 441 native PORT
dispersion fits consumed 1.496 s cumulative. The nested objective and
intercept entries cannot be added together. The evidence still identifies
parameter fitting—not residual generation—as the first optimization target
for this fixture.

Machine-readable evidence is in `results/benchmark.json`,
`results/quality.json`, `results/python_profile_top.json` and
`results/environment.json`.
