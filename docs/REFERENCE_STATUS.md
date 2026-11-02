# v0.4.0 reference status

## Purpose

Version 0.4.0 completes the preliminary SCTransform v2 reference work in this
repository. The package remains installable and useful as a CPU
implementation, while its strict oracle and profiling evidence provide the
baseline for a separate MIT-compatible native Python CPU/GPU project.

## Pinned statistical contract

| component | version / source |
|---|---|
| R | 4.5.2, `b4322f3207e8d99d9dc6e6fa6b80ad23c4dc2f5f` |
| sctransform | 0.4.3, `49e35b5aeb76a602910207cbfde1561093340be3` |
| glmGamPoi | 1.22.0, `908cde9d595259b4a317343ddfbcf20036ef0895` |
| matrixStats | 1.5.0 |

Validated behavior includes the default fixed-slope v2 path, R-compatible
sampling, numeric regularized covariates, one categorical batch variable,
finite-theta nonreg covariates, Pearson residuals, residual variance, variable
feature ordering and corrected UMI for the supported designs.

Alternative fit methods and residual types, arbitrary R contrasts, complete
Seurat object semantics, macOS and unvalidated platform math libraries are
outside the contract. Official nonreg fitting also fails when regularization
produces `theta=Inf`; Python fails closed for the same condition.

## Native artifacts

| artifact | build contract | SHA-256 |
|---|---|---|
| `sct_port_r45.dll` | Rtools45 GCC/GFortran 14.3.0 | `e2abaa8f7d22ec193e48c9b75f3294e56f763e7a4f7e2318150a292ad5b10ba2` |
| `libsct_port_r45.so` | Ubuntu 26.04, GCC/GFortran 15.2.0, glibc 2.43 | `97b426f7beced380fd98ded949cd1ad736decff503638d29777616efaf1b3645` |

The source distribution contains the corresponding source and build scripts.
R is used only by the development oracle and is not a runtime dependency.

## Correctness evidence

- the 800 x 240 PBMC3K-derived fixture is binary-identical to fresh,
  same-platform R oracles on the two pinned environments;
- the adversarial fixture covers low counts, the Poisson boundary, clipping,
  rare-gene filtering and an extreme library;
- numeric, batch, corrected and batch-plus-nonreg live gates compare complete
  parameters, residuals, residual statistics and corrected integers;
- functional CI covers Python 3.10--3.13 on Windows and Linux;
- Windows and Linux wheels are installed in temporary environments and load
  their packaged native libraries.

## Performance evidence

The controlled small benchmark used one warm-up and five fresh-process runs
with one BLAS/OpenMP thread. On this fixture:

| platform | Python compute | R compute | ratio | Python/R sampled peak RSS |
|---|---:|---:|---:|---:|
| Windows | 0.491 s | 3.380 s | 6.88x | 63.2 / 731.3 MiB |
| Linux/WSL2 | 0.206 s | 3.185 s | 15.46x | 64.7 / 680.0 MiB |

These measurements demonstrate a real acceleration on the named fixture and
machines. They are not a claim about larger data, another platform or every
SCT configuration. The profile identified per-gene parameter fitting as the
dominant small-fixture CPU stage; broad scaling remains work for the separate
native implementation.

## Handoff

Future work on a dependency-clean NumPy/SciPy CPU implementation, AnnData,
CuPy/RMM, streaming, async execution and GPU kernels moves to a separate
MIT-compatible repository. This GPL package remains the external reference
and differential-testing oracle; its R/Fortran-derived code is not an
implementation dependency of that project.
