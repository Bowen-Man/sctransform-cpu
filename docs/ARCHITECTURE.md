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
per-batch SCT v2 sequence. `r_glmgampoi.py`, `r45_mathlib.py` and the bundled
PORT library form the numerical compatibility boundary. The release has no
GPU backend, asynchronous execution or experimental solver selection.

Residuals are generally dense even when input UMI counts are sparse because a
zero count usually has a nonzero Pearson residual. The implementation keeps
the input CSR and densifies one gene tile at a time.
