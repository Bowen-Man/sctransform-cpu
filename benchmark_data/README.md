# Benchmark fixture provenance

This is an 800-gene x 240-cell computational subset of the raw UMI PBMC3K
matrix used by the Seurat SCTransform vignette. The source workflow is
documented in the official
[Seurat SCTransform vignette](https://github.com/satijalab/seurat/blob/main/vignettes/sctransform_vignette.Rmd).

Construction retained the first 240 cells and the first 800 genes detected in
at least five selected cells. Gene identifiers in this public regression
fixture are deliberately replaced by stable `gene-xxxxx` titles; it must not
be used for biological interpretation. The matrix has 22,131 nonzero entries.

`counts.mtx` and `counts.npz` encode the same gene x cell integer matrix for R
and Python respectively. `genes.tsv`/`cells.tsv` and the corresponding text
files preserve identical row and column order.
