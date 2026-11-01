# Code map and development rules

## Runtime modules

| Module | Responsibility | Important invariant |
|---|---|---|
| `core.py` | Public API, filtering, stage order, tiling and result construction | Input is raw gene x cell integer UMI; no silent subsampling |
| `fit.py` | Ordered batches of at most 500 genes | Preserve R batch and gene order |
| `model.py` | Moment start and one-intercept refits | Preserve scalar accumulation and convergence order |
| `r_glmgampoi.py` | Dispersion objective, score, Hessian and fallback chain | Match pinned glmGamPoi 1.22 semantics |
| `port_native.py` | Typed `ctypes` boundary for the bundled optimizer | Fail closed if the native library is missing or returns an invalid state |
| `r45_mathlib.py` | Audited R 4.5-compatible scalar functions | No silent SciPy/libm substitution in the strict path |
| `r_compat_stats.py` | Sparse row statistics, bandwidth, outliers and kernel smoothing | Cross-gene order is part of the numerical contract |
| `residuals.py` | Float64 Pearson residual tiles and sample moments | Sparse zeros have nonzero residuals; output is generally dense |
| `corrected_counts.py` | Reverse unclipped residuals at the target depth | Integer rounding order is part of the contract |
| `native/` | R PORT Fortran, compatibility callbacks and Windows DLL | Corresponding source and provenance must ship with releases |

The public package deliberately contains no GPU backend, async executor,
alternative solver switch or hidden approximate fallback.

## Research and release layers

- `benchmarks/bit_reference.R` writes raw little-endian IEEE-754 values; it
  avoids CSV formatting and is the authoritative bit-level oracle.
- `benchmarks/bit_identity.py` compares float bit patterns, integers and
  ordered titles. `--enforce` exits nonzero on one mismatched bit.
- `benchmarks/compare.py` is the less strict numerical acceptance gate.
- `benchmarks/run_benchmark.py` measures identical work in fresh Python and R
  subprocesses, including sampled RSS and startup separately from compute.
- `tests/data/solver_stress*` is the stored 128-case dispersion regression.
- `results/` stores machine-readable evidence; prose reports summarize rather
  than replace it.

## Development principles

1. **Pin the contract.** Every claim names the R, sctransform, glmGamPoi,
   operating-system and compiler versions it actually covers.
2. **Fail closed.** Missing selection manifests, native code or unsupported
   model features raise an error instead of choosing a nearby algorithm.
3. **Use layered gates.** Bit identity, numerical tolerance, statistical
   ranking and performance are different questions and are reported
   separately.
4. **Compare equal work.** Counts, ordering, fitted cells/genes, returned
   residuals, corrected counts and thread settings must match.
5. **Keep R out of the runtime.** R is an oracle and benchmark dependency, not
   an end-user dependency.
6. **Preserve numerical order before optimizing.** Reassociation, vectorized
   reductions, sparse-format changes and fusion require new oracle evidence.
7. **Measure end to end.** A faster inner callback is not an SCT speedup until
   stage and process wall time, RSS and output quality are measured.
8. **Optimize behind the CPU reference.** GPU or alternative CPU backends may
   be added only after the strict CPU contract is stable.
9. **Keep provenance releasable.** GPL notices, corresponding native source,
   fixture provenance, build commands and exact hashes are release artifacts.

## Current bit-identity status

The first binary oracle run on the checked-in 800 x 240 PBMC fixture is red.
The earlier tolerance gate remains green, but it must not be described as bit
identity.

| Value | elements | bit mismatches |
|---|---:|---:|
| library size | 240 | 0 |
| arithmetic mean | 800 | 0 |
| sample variance | 800 | 0 |
| log geometric mean | 800 | 1 |
| step-1 intercept | 441 | 3 |
| step-1 theta | 441 | 40 |
| step-1 fixed slope | 441 | 0 |
| regularized intercept | 800 | 400 |
| regularized theta | 800 | 439 |
| regularized slope | 800 | 387 |
| residual mean | 800 | 800 |
| residual variance | 800 | 749 |
| Pearson residual | 192,000 | 177,065 |
| outlier mask | 441 | 0 |
| corrected UMI | 192,000 | 0 |

Gene/cell/fit ordering is identical. The complete residual-variance order is
not bit-identical even though the previously reported Top-100/200/300 sets are
identical.

Run the hard gate locally:

```powershell
$env:SCT_TEST_RSCRIPT = 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe'
python -m pytest -q tests/test_r_bit_identity.py
```

The test is intentionally an ordinary assertion, not an `xfail`. It currently
fails and will become a required release gate only after every compared field
passes.

## Remaining work, in order

### P0: pinned Windows bit identity

1. Make the final log-geometric-mean case use the exact R exponent/log path.
2. Move the remaining glmGamPoi start, intercept and callback arithmetic into
   one audited native boundary; 40 theta trajectories still diverge.
3. Port the exact R/sctransform bandwidth, outlier and `ksmooth` operation
   order, including the nominally constant slope smoothing.
4. Match the residual formula and matrixStats reduction order after the model
   parameters are bit-identical.
5. Add more binary oracles before treating one PBMC fixture as a guarantee.

This guarantee can only be version/platform pinned. Cross-platform bit
identity is not a reasonable default contract because libm, compiler and
extended-precision behavior differ; cross-platform releases should retain a
separate numerical-tolerance contract.

### P1: statistical feature coverage

- exact independent R RNG/density-weighted sampling, including ties;
- batch and additional covariates;
- complete Seurat-facing output semantics;
- more real datasets, extreme libraries, Poisson boundaries and rare counts.

### P2: packaging and scale

- Python 3.10, 3.11 and 3.13 Windows wheels;
- Linux/macOS native builds with their own numerical gates;
- PBMC10K/30K and controlled larger scaling, peak allocation timelines and
  thread-scaling measurements;
- AnnData/Scanpy integration without changing the core matrix contract.

### P3: optimization backends

Only after P0/P1 gates are explicit: vectorized/native CPU callbacks,
streamed residual reductions, then evidence-driven GPU prototypes. Current
profiling identifies per-gene fitting as the first CPU hotspot; it does not
yet justify moving the whole pipeline to GPU.
