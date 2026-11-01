# SCTransform v2 source map and compatibility contract

## Scope and pinned sources

This document maps the direct `sctransform::vst(..., vst.flavor="v2")`
computation. It does not claim compatibility with Seurat object construction,
assay replacement, integration preparation, or every historical SCT method.

The audit is pinned to:

| component | version / commit | role |
|---|---|---|
| sctransform | 0.4.3, `49e35b5aeb76a602910207cbfde1561093340be3` | authoritative workflow |
| R | 4.5.2; source branch snapshot `b4322f3207e8d99d9dc6e6fa6b80ad23c4dc2f5f` | RNG, `density`, `sample`, model-matrix semantics |
| glmGamPoi | 1.22.0; source snapshot `908cde9d595259b4a317343ddfbcf20036ef0895` | Gamma-Poisson fitting used by v2 |
| matrixStats | 1.5.0 | row statistics used by the oracle |

Line references below refer to the pinned
[`R/vst.R`](https://github.com/satijalab/sctransform/blob/49e35b5aeb76a602910207cbfde1561093340be3/R/vst.R)
and
[`R/fit.R`](https://github.com/satijalab/sctransform/blob/49e35b5aeb76a602910207cbfde1561093340be3/R/fit.R).
Statements labelled **confirmed** come directly from those sources. Entries
labelled **implemented** have a Python gate. **Missing** means the public
Python API must reject or omit the branch rather than silently approximate it.

## Flavor rewrite

`vst` first receives the general defaults at `vst.R:109-138`. The v2 flavor
then rewrites four decisions at `vst.R:139-159`:

1. choose `glmGamPoi_offset` when glmGamPoi is installed, otherwise
   `nb_offset`;
2. force `exclude_poisson=TRUE`;
3. change `min_variance=-Inf` to `"umi_median"`;
4. change a null `n_cells` to 2,000.

It does **not** rewrite `n_genes`, whose general default remains 2,000. The
separate `startsWith(method, "offset")` rewrite at `vst.R:173-187` does not
match the string `glmGamPoi_offset`; treating those as the same branch is a
source-reading error.

## End-to-end data flow

```text
gene x cell raw UMI (sparse or dense)
  -> make/validate cell attributes
  -> retain genes detected in at least min_cells cells
  -> full-data log geometric mean
  -> optional uniform cell sample
  -> recheck min_cells on sampled cells
  -> exclude full-data variance <= mean genes from step 1
  -> optional inverse-density weighted gene sample
  -> construct model formula/design
  -> fit per-gene Gamma-Poisson parameters in bins
  -> clamp theta >= 1e-7
  -> detect parameter outliers and regularize across genes
  -> optional non-regularized covariate fit
  -> generate dense Pearson/deviance residuals in gene bins
  -> optional corrected counts
  -> clip residuals
  -> optional gene attributes and return metadata
```

Source anchors are `vst.R:192-292` for preprocessing and fitting,
`vst.R:294-345` for regularization/non-regularized variables,
`vst.R:348-407` for residuals, and `vst.R:409-463` for construction of the
returned object.

## Parameter and branch matrix

| official argument | v2 behavior confirmed from source | Python contract |
|---|---|---|
| `umi` | genes x cells; gene and cell names are operational identifiers | **implemented**, CSR canonicalization, raw integer validation |
| `cell_attr` | created/validated by `make_cell_attr`; may supply arbitrary requested covariates | **missing** for arbitrary covariates |
| `latent_var` | default `log_umi`; additional variables enter the design | only `log_umi` **implemented** |
| `batch_var` | factorized; formula becomes `(latent):batch + batch + 0` | **missing**, independently gated work item |
| `latent_var_nonreg` | fitted after regularization and appended to final design | **missing** |
| `n_genes` | default 2,000; inverse-density sample when fewer than eligible genes | R 4.5 RNG and density sampler **implemented**, live-R selection gate required |
| `n_cells` | null becomes 2,000 under v2; uniform sample without replacement | R 4.5 RNG sampler **implemented** |
| `method` | v2 selects `glmGamPoi_offset`, fallback `nb_offset` | pinned glmGamPoi offset path **implemented**; fallback **missing** |
| `do_regularize` | default true | true **implemented**; false **missing** |
| `theta_regularization` | default `od_factor`; `log_theta` is an alternative | `od_factor` **implemented** |
| `res_clip_range` | default `[-sqrt(C), sqrt(C)]`, applied after residual construction | default **implemented**; custom bounds **missing** |
| `bin_size` | controls fit and residual gene bins | fit 500 and configurable residual tile **implemented**; exact general API missing |
| `min_cells` | detection filter before and after cell sampling | **implemented** |
| `residual_type` | `pearson`, `deviance`, or `none` | Pearson **implemented**; other branches **missing** |
| `return_cell_attr` | removes `cell_attr` when false | compact result has no cell attributes |
| `return_gene_attr` | detection rate, gmean, amean, variance and residual summaries | **implemented** in `SCTResult.to_vst_dict()` |
| `return_corrected_umi` | only with Pearson residuals; rounded positive sparse counts | **implemented** for targeted path |
| `min_variance` | v2 default `umi_median`; also numeric, `model_mean`, `model_median` | `umi_median` **implemented** |
| `bw_adjust` | scales `bw.SJ` during cross-gene smoothing | **implemented** |
| `gmean_eps` | default 1 in sparse geometric mean | fixed value 1 **implemented** |
| `theta_estimation_fun` | used by non-glmGamPoi methods | outside current glmGamPoi contract |
| `theta_given` | used by given/shared-theta methods | outside current contract |
| `exclude_poisson` | forced true by v2 | **implemented** for the unbatched path |
| `use_geometric_mean` | controls abundance coordinate for sampling/smoothing | true **implemented** |
| `use_geometric_mean_offset` | affects fixed-parameter branches | false only |
| `fix_intercept`, `fix_slope` | affect non-offset/general fitting paths | not exposed; v2 offset slope is fixed by its own fitting function |
| `scale_factor` | corrected-count target depth; NA uses median UMI | default and numeric scale factor **implemented** |
| `vst.flavor` | `"v2"` performs the rewrite above | fixed to v2 |
| `verbosity` | messages/progress only | deliberately omitted from numerical core |

## Filtering and sampling semantics

### Gene and cell filtering

- Initial detection is `rowSums(umi >= 0.01) >= min_cells`
  (`vst.R:198-202`). For integer UMI this is equivalent to nonzero count.
- After a cell sample, a gene must again be positive in at least `min_cells`
  sampled cells (`vst.R:216-230`).
- v2 Poisson exclusion uses the **full retained matrix** mean and sample
  variance, keeping `variance - mean > 0` for step 1
  (`vst.R:241-258`).

### Cell sample

`sample(colnames(umi), n_cells)` is uniform without replacement and preserves
the random draw order. Python reproduces R 4.5's Mersenne-Twister seed
scramble, rejection integer sampler, and swap-with-tail algorithm. The seed is
an explicit Python configuration value because R's seed lives outside `vst`.

### Gene sample

At `vst.R:262-273`, official code:

1. estimates `density(log_gmean, bw="nrd")` on a 512-point grid;
2. linearly interpolates density back to each gene;
3. assigns probability `1 / (density + double.eps)`;
4. samples without replacement using R's weighted sampler.

The Python path ports the R 4.5 binned-density/FFT data flow, `revsort`
heapsort and sequential weighted sampler. The contract is the ordered selected
gene titles, not binary identity of an internal FFT grid.

## v2 expression model

For the unbatched default:

\[
Y_{gc}\sim\mathrm{NB}(\mu_{gc},\theta_g),\qquad
\operatorname{Var}(Y_{gc})=\mu_{gc}+\mu_{gc}^2/\theta_g,
\]

\[
\log\mu_{gc}=\beta_{0g}+\log U_c
=\beta_{0g}+\log(10)\log_{10}U_c.
\]

This parameterization is **confirmed** by `fit.R:107-176`. The official
function creates natural-log library size as an offset, removes `log_umi` from
the no-batch formula, fits the remaining coefficients with `glm_gp`, and then
re-inserts a `log(10)` coefficient under the name `log_umi`.

With a batch variable, the source detects `(log_umi) : batch`, sets the offset
to zero, and lets glmGamPoi estimate the batch-specific columns. The formula
is assembled at `vst.R:275-279`. Therefore the batch path is not the same
fixed-slope solver repeated by batch.

## glmGamPoi fitting sequence

For each gene bin, `get_model_pars` converts the sparse slice to a dense matrix
and calls `fit_glmGamPoi_offset` (`vst.R:551-630`). glmGamPoi 1.22 then:

1. creates the offset matrix and model matrix;
2. obtains rough dispersion and coefficient estimates;
3. fits coefficients by group-wise Newton/Fisher scoring when the design has
   a small number of distinct rows, otherwise general Fisher scoring;
4. maximizes Gamma-Poisson dispersion with Cox-Reid adjustment;
5. computes the dispersion trend used for its internal shrinkage;
6. refits coefficients using the latest dispersion.

Genes are independent during coefficient/dispersion fitting, but share the
bin-level dispersion trend and later SCT regularization. This is why “each
gene is independent” is only conditionally true.

## Poisson boundary and SCT regularization

Official v2 excludes genes with full-data `variance <= mean` from step-1
fitting. During regularization it also defines Poisson genes as
`variance <= mean` or `mean < 1e-3`; their final theta is infinite and fixed
mean-model coefficients are used. Fitted genes with extremely discrepant
naive and MLE theta are excluded from smoothing. The current unbatched Python
regularizer reproduces this pinned path.

The batch path is materially different: `reg_model_pars` disables the same
fixed-slope/Poisson shortcut for batch-specific coefficients. It requires a
coefficient-wise smoothing and design-matrix oracle before it can be marked
implemented.

## Residual and corrected-count contract

For Pearson output:

\[
r_{gc}=\frac{Y_{gc}-\hat\mu_{gc}}
{\sqrt{\max(\hat\mu_{gc}+\hat\mu_{gc}^2/\hat\theta_g, v_{min})}}.
\]

For the v2 default, `v_min=(median(nonzero UMI)/5)^2`
(`vst.R:348-360`). Residuals are dense even when UMI is sparse because a zero
count normally has a nonzero negative residual. The official main `vst` path
allocates the full dense gene x cell result (`vst.R:370`); its helper
`get_residual_var` can instead reduce residual variance by bins without
retaining the full matrix.

The Python `output="hvg"` path similarly computes tiled statistics and only
materializes selected rows. `output="full"` is required for the R-like return
adapter. Corrected counts are reconstructed from Pearson residuals, clipped
to nonnegative values, rounded, and returned sparse.

## Return schema

The official list is constructed at `vst.R:409-419`, then conditionally gains
`umi_corrected`, `gene_attr`, and `times` at `vst.R:423-463`. The stable names
are:

```text
y, model_str, model_pars, model_pars_outliers, model_pars_fit,
model_str_nonreg, model_pars_nonreg, arguments,
genes_log_gmean_step1, cells_step1, cell_attr,
[umi_corrected], [gene_attr], times
```

`SCTResult.to_vst_dict()` exposes these names using Python containers and
explicit row/column title arrays. It is a structural adapter, not an emulation
of an R list or Seurat assay.

## Alignment gates and unresolved work

| gate | exact target | status |
|---|---|---|
| pinned unbatched no-sampling fixture | binary64/int32/title identity | green |
| R 4.5 uniform RNG and ordered sample | exact values/titles | green |
| inverse-density gene selection | ordered titles against live R | new gate |
| Poisson equality, low-mean, clipping, rare genes, zero filtering | branch/output identity | new adversarial suite |
| return names, shapes, titles and optional fields | structural identity | new contract tests |
| numeric covariate | model matrix, coefficients, residuals | missing |
| categorical batch | factor levels, coefficients, regularization, residuals | missing |
| non-regularized variables | appended design and coefficients | missing |
| fallback `nb_offset` | parameters and failure behavior | missing |

No missing row may be advertised as supported. In particular, a tolerant fit
from statsmodels is not evidence of glmGamPoi compatibility, and a matching
HVG list is not evidence of parameter or residual identity.
