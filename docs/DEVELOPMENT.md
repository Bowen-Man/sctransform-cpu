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
- `benchmarks/profile_python.py` performs one warm-up and then profiles only
  `sctransform`, keeping import, I/O and DLL initialization outside cProfile.
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

The checked-in 800 x 240 PBMC fixture is green under the pinned Windows
contract. Every compared value has zero mismatches:

- library size, offsets and `log_umi`;
- arithmetic/geometric means and sample variance;
- glmGamPoi starts, fitted intercept, theta and fixed slope for 441 genes;
- bandwidth, outlier mask and all regularized parameters;
- all 192,000 model means, model variances, residual numerators,
  denominators and Pearson residuals;
- all residual means and variances, complete HVG order and 192,000 corrected
  UMI integers;
- retained genes/cells and ordered fit-gene/fit-cell manifests.

The machine-readable evidence is `results/bit_identity.json`. A passing live
test also regenerated the oracle with R 4.5.2 rather than reusing that file.

Run the hard gate locally:

```powershell
$env:SCT_TEST_RSCRIPT = 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe'
python -m pytest -q tests/test_r_bit_identity.py
```

The test is intentionally an ordinary assertion, not an `xfail`. It is the
required Windows release gate and fails on the first mismatched bit or title.

## Remaining work, in order

### P0: harden the pinned Windows guarantee

1. Add independent binary oracles covering low-count, Poisson-boundary,
   clipping and extreme-library cases; one PBMC fixture is a regression gate,
   not a universal proof.
2. Build and test Windows wheels for Python 3.10, 3.11, 3.12 and 3.13.
3. Record compiler, DLL and oracle hashes as release metadata and reproduce
   the native build in CI.
4. Exercise different tile sizes and requested output modes against the same
   full R result.

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
