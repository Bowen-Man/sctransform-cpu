args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop('usage: Rscript selection_reference.R DATA OUTPUT SEED N_CELLS N_GENES')
}
suppressPackageStartupMessages({
  library(Matrix)
  library(sctransform)
  library(jsonlite)
})
expected <- c(R='4.5.2', sctransform='0.4.3', glmGamPoi='1.22.0')
actual <- c(
  R=as.character(getRversion()),
  sctransform=as.character(packageVersion('sctransform')),
  glmGamPoi=as.character(packageVersion('glmGamPoi'))
)
if (!identical(actual, expected)) stop('strict selection oracle version mismatch')

data_dir <- args[1]
counts <- as(readMM(file.path(data_dir, 'counts.mtx')), 'CsparseMatrix')
rownames(counts) <- read.delim(file.path(data_dir, 'genes.tsv'))$gene
colnames(counts) <- read.delim(file.path(data_dir, 'cells.tsv'))$cell
set.seed(as.integer(args[3]))
result <- vst(
  counts,
  vst.flavor='v2',
  n_cells=as.integer(args[4]),
  n_genes=as.integer(args[5]),
  min_cells=5,
  residual_type='none',
  return_gene_attr=FALSE,
  verbosity=0
)
write_json(
  list(cells=unname(result$cells_step1), genes=unname(rownames(result$model_pars))),
  args[2],
  auto_unbox=TRUE,
  pretty=TRUE
)
