args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('usage: Rscript bit_reference.R DATA OUTPUT_DIR')
suppressPackageStartupMessages({library(Matrix); library(sctransform); library(jsonlite)})
data_dir <- args[1]
output <- args[2]
if (dir.exists(output) && length(list.files(output, all.files=TRUE, no..=TRUE))) {
  stop('output directory is not empty')
}
dir.create(output, recursive=TRUE, showWarnings=FALSE)

counts <- as(readMM(file.path(data_dir, 'counts.mtx')), 'CsparseMatrix')
rownames(counts) <- read.delim(file.path(data_dir, 'genes.tsv'))$gene
colnames(counts) <- read.delim(file.path(data_dir, 'cells.tsv'))$cell
set.seed(42L)
result <- vst(
  counts,
  vst.flavor='v2',
  n_cells=ncol(counts),
  n_genes=nrow(counts),
  min_cells=5,
  return_gene_attr=TRUE,
  return_corrected_umi=TRUE,
  verbosity=0
)

write_f64 <- function(name, value, matrix=FALSE) {
  connection <- file(file.path(output, paste0(name, '.f64le')), 'wb')
  on.exit(close(connection), add=TRUE)
  payload <- if (matrix) as.double(t(value)) else as.double(value)
  writeBin(payload, connection, size=8L, endian='little')
}
write_i32 <- function(name, value, matrix=FALSE) {
  connection <- file(file.path(output, paste0(name, '.i32le')), 'wb')
  on.exit(close(connection), add=TRUE)
  payload <- if (matrix) as.integer(t(value)) else as.integer(value)
  writeBin(payload, connection, size=4L, endian='little')
}

writeLines(rownames(result$model_pars), file.path(output, 'fit_genes.txt'))
writeLines(result$cells_step1, file.path(output, 'fit_cells.txt'))
writeLines(rownames(result$model_pars_fit), file.path(output, 'genes.txt'))
writeLines(colnames(result$y), file.path(output, 'cells.txt'))
writeLines(
  rownames(result$model_pars_fit)[order(-result$gene_attr$residual_variance, method='radix')],
  file.path(output, 'hvg_order.txt')
)

write_f64('library', Matrix::colSums(counts))
write_f64('amean', result$gene_attr$amean)
write_f64('variance', result$gene_attr$variance)
write_f64('log_gmean', log10(result$gene_attr$gmean))
write_f64('fit_intercept', result$model_pars[, '(Intercept)'])
write_f64('fit_theta', result$model_pars[, 'theta'])
write_f64('fit_slope', result$model_pars[, 'log_umi'])
write_i32('fit_outlier', result$model_pars_outliers)
write_f64('intercept', result$model_pars_fit[, '(Intercept)'])
write_f64('theta', result$model_pars_fit[, 'theta'])
write_f64('slope', result$model_pars_fit[, 'log_umi'])
write_f64('residual_mean', result$gene_attr$residual_mean)
write_f64('residual_variance', result$gene_attr$residual_variance)
write_f64('residual', result$y, matrix=TRUE)
write_i32('corrected_umi', as.matrix(result$umi_corrected), matrix=TRUE)

metadata <- list(
  byte_order='little',
  float_format='IEEE-754 binary64',
  retained_genes=nrow(result$model_pars_fit),
  fitted_genes=nrow(result$model_pars),
  fitted_cells=length(result$cells_step1),
  cells=ncol(result$y),
  r_version=as.character(getRversion()),
  sctransform_version=as.character(packageVersion('sctransform')),
  glmgampoi_version=as.character(packageVersion('glmGamPoi'))
)
write_json(metadata, file.path(output, 'metadata.json'), auto_unbox=TRUE, pretty=TRUE)
