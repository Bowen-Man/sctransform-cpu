# Architecture

```text
raw gene x cell UMI CSR
  -> library size before filtering
  -> gene filtering and R-order row statistics
  -> explicit fit-cell / fit-gene selection
  -> batches of at most 500 genes
       -> moment dispersion and intercept
       -> fixed-mean glmGamPoi dispersion MLE
       -> local dispersion trend
       -> intercept refit
  -> cross-gene outlier detection and kernel regularization
  -> tiled float64 Pearson residual statistics
  -> residual-variance ranking
  -> optional residual or corrected-count materialization
```

`core.py` owns the data flow and public contract. `fit.py` contains the
per-batch SCT v2 sequence. `r_glmgampoi.py`, `r45_mathlib.py`,
`r_compat_stats.py`, `port_native.py` and the bundled native library form the
numerical compatibility boundary. The library contains the PORT optimizer
and narrowly ported R/Rcpp, glmGamPoi and matrixStats arithmetic; it does not
embed or call an R runtime. The release has no GPU backend, asynchronous
execution or experimental solver selection.

Residuals are generally dense even when input UMI counts are sparse because a
zero count usually has a nonzero Pearson residual. The implementation keeps
the input CSR and densifies one gene tile at a time.

The exact path preserves numerical operation order as part of the API
contract. In particular, residual means reproduce `base::rowMeans` long-double
accumulation, while residual variances reproduce matrixStats 1.5.0's binary64
refinement pass. Replacing either with a NumPy reduction changes bits even
when the statistical result is numerically close.
