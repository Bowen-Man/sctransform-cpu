# Python and RAPIDS roadmap

## Mission

Build an R-free Python implementation of the SCTransform v2 computation, use
the CPU implementation as the statistical and numerical reference, optimize
that reference from profiler evidence, and only then design an optional GPU
backend for AnnData/CuPy/RAPIDS-oriented workflows.

"R-free" describes the installed runtime. R remains the authoritative oracle
used by development and release tests. The existing small native compatibility
library is also part of the Python package: it reproduces platform arithmetic
that NumPy/SciPy alone do not match bit-for-bit.

The primary parity target is the SCTransform v2 computation exposed by
`sctransform::vst`. Reproducing the complete Seurat R object lifecycle is a
separate integration concern and must not be confused with statistical parity.

## Current compatibility state

| Area | State | Evidence or gap |
|---|---|---|
| R-free Python runtime | complete | R is absent from package dependencies and runtime calls |
| Direct unbatched v2 path with fixed `log_umi` slope | validated | live Windows and Linux R 4.5.2 binary oracles |
| Model parameters, Pearson residuals and residual variance | validated in current contract | zero mismatches on the pinned 800 x 240 fixture |
| Corrected UMI counts | validated in current contract | integer-exact comparison on the pinned fixture |
| Explicit shared fit-cell and fit-gene selections | validated | ordered manifests are compared with R |
| R-compatible implicit sampling | validated in pinned contract | live R selection gates compare complete ordered cell/gene titles |
| Regularized numeric covariates | validated | live Windows and Linux R binary oracles compare coefficients and residual outputs |
| One categorical batch variable | validated | live Windows and Linux R binary oracles cover batch-specific intercepts/slopes and regularization |
| General-design corrected UMI | validated | numeric-covariate and batch+nonreg live binary oracles on Windows and Linux |
| Non-regularized latent variables | validated with finite theta | batch+nonreg live oracle; the official theta=Inf failure is reproduced by failing closed |
| Complete v2 flavor/method option coverage | incomplete | source map exists; alternative method/residual branches remain outside the contract |
| Poisson-boundary, clipping, rare-count and extreme-library suites | validated for current path | adversarial live-R fixture is bit-identical; broader property coverage remains ongoing |
| Targeted `vst` return schema | implemented | `to_vst_dict()` exposes official names with explicit matrix titles; unsupported branches remain absent |
| Broad real-data and scaling evidence | incomplete | current benchmark is intentionally small |
| AnnData/CuPy/RAPIDS interface | not implemented | deferred until CPU gates are complete |
| GPU backend | not implemented | no GPU speedup claim is currently made |

## Phase 1: complete the alignment contract

### Work

1. Pin the official R source commit and package/runtime matrix used by every
   oracle.
2. Produce a source-level option map for every targeted v2 branch: defaults,
   filters, sampling, model fitting, Poisson handling, regularization,
   residual selection, clipping, corrected counts and returned metadata.
3. Implement missing branches one at a time behind fixtures that share counts,
   names, cell/gene selections and complete parameter manifests with R.
4. Add independent real and adversarial datasets for low counts, Poisson
   boundaries, clipping, extreme library sizes, rare genes and all-zero
   filtering behavior.
5. Separate exact binary gates from documented numerical/statistical gates.
   A feature is not aligned merely because its HVG list looks similar.
6. Define a stable Python result schema and an optional AnnData adapter without
   forcing AnnData into the numerical core.

### Exit gate

- Every targeted `vst` v2 option is either green against a live pinned R oracle
  or explicitly declared outside the package contract.
- No unknown or silently approximated branch remains in the public API.
- Windows and Linux release fixtures pass their own platform contracts.
- The compatibility table and public documentation match the implementation.

## Phase 2: CPU optimization and implementation hardening

### Work

1. Freeze the Phase 1 oracle suite before changing execution order.
2. Benchmark small, medium and resource-bounded large workloads across cells,
   genes, density, output mode and thread configuration.
3. Record end-to-end and stage wall time, CPU time, peak RSS, native
   allocations, matrix conversions and intermediate lifetimes.
4. Classify measured hotspots as solver, Python dispatch, allocation, sparse
   traversal, vector math, memory bandwidth or algorithmic complexity.
5. Evaluate low-risk dataflow/allocation changes before compiled or parallel
   alternatives. Every optimized path retains the CPU reference as an oracle.
6. Report kernel/function, stage and end-to-end effects separately; reject
   changes whose local gains do not survive end-to-end measurement.

### Exit gate

- Reproducible baselines and profiles exist for the agreed workload matrix.
- Each accepted optimization has correctness, time and memory evidence.
- CPU reference and optimized modes share the same public semantics.
- Remaining bottlenecks and scaling laws are quantified rather than inferred.

## Phase 3: GPU design after CPU evidence

Only after both previous gates pass will the project publish a GPU feasibility
and implementation proposal. That proposal must:

- identify the exact stages worth accelerating and the end-to-end upper bound;
- compare NumPy/SciPy, Numba, CuPy, cuSPARSE, Triton, CUDA C++ and CPU/GPU
  hybrid execution without selecting a winner in advance;
- account for H2D/D2H transfers, allocation, synchronization and conversion
  between SciPy, CuPy and AnnData representations;
- decide sparse, dense-tiled and streaming behavior per stage rather than for
  the entire pipeline;
- estimate working-set and VRAM requirements for a 12 GB consumer GPU;
- preserve the Phase 1 oracle suite for every backend;
- define a narrow backend protocol so the CPU implementation remains usable
  without RAPIDS or CUDA dependencies.

The first GPU deliverable is therefore an evidence-backed design and isolated
microbenchmark set, not a premature all-GPU rewrite.
