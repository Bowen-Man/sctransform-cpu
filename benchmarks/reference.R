args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('usage: Rscript reference.R DATA OUTPUT_DIR')
suppressPackageStartupMessages({library(Matrix); library(sctransform)})
script_path <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_path)), 'r_write_exact_csv.R'))
data_dir <- args[1]
output <- args[2]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
if (length(list.files(output, all.files=TRUE, no..=TRUE))) stop('output is not empty')
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
options(digits=17)
write_exact_csv(
  data.frame(gene=rownames(result$model_pars), result$model_pars,
             outlier=result$model_pars_outliers, check.names=FALSE),
  file.path(output, 'fit.csv'), row.names=FALSE)
write_exact_csv(
  data.frame(gene=rownames(result$model_pars_fit), result$model_pars_fit,
             result$gene_attr, check.names=FALSE),
  file.path(output, 'genes.csv'), row.names=FALSE)
write_exact_csv(
  data.frame(gene=rownames(result$y), result$y, check.names=FALSE),
  file.path(output, 'residuals.csv'), row.names=FALSE)
write_exact_csv(
  data.frame(gene=rownames(result$umi_corrected), as.matrix(result$umi_corrected),
             check.names=FALSE),
  file.path(output, 'corrected_umi.csv'), row.names=FALSE)
writeLines(result$cells_step1, file.path(output, 'fit_cells.txt'))
capture.output(sessionInfo(), file=file.path(output, 'session_info.txt'))
