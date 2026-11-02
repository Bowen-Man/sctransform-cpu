args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('usage: Rscript design_reference.R DATA OUTPUT MODE')
suppressPackageStartupMessages({library(Matrix); library(sctransform)})
expected <- c(R='4.5.2', sctransform='0.4.3', glmGamPoi='1.22.0', matrixStats='1.5.0')
actual <- c(R=as.character(getRversion()),
            sctransform=as.character(packageVersion('sctransform')),
            glmGamPoi=as.character(packageVersion('glmGamPoi')),
            matrixStats=as.character(packageVersion('matrixStats')))
if (!identical(actual, expected)) stop('strict oracle version mismatch')

counts <- as(readMM(file.path(args[1], 'counts.mtx')), 'CsparseMatrix')
rownames(counts) <- read.delim(file.path(args[1], 'genes.tsv'))$gene
colnames(counts) <- read.delim(file.path(args[1], 'cells.tsv'))$cell
cell_attr <- data.frame(
  percent=((seq_len(ncol(counts))-1) %% 17) / 17,
  batch=rep(c('A', 'B'), length.out=ncol(counts)),
  row.names=colnames(counts)
)
result <- switch(args[3],
  covariate=vst(counts, vst.flavor='v2', latent_var=c('log_umi', 'percent'),
                cell_attr=cell_attr, n_cells=ncol(counts), n_genes=nrow(counts),
                min_cells=5, return_gene_attr=TRUE, verbosity=0),
  covariate_corrected=vst(
                counts, vst.flavor='v2', latent_var=c('log_umi', 'percent'),
                cell_attr=cell_attr, n_cells=ncol(counts), n_genes=nrow(counts),
                min_cells=5, return_gene_attr=TRUE, return_corrected_umi=TRUE,
                scale_factor=1234, verbosity=0),
  batch=vst(counts, vst.flavor='v2', batch_var='batch', cell_attr=cell_attr,
            n_cells=ncol(counts), n_genes=nrow(counts), min_cells=5,
            return_gene_attr=TRUE, verbosity=0),
  batch_nonreg=vst(
            counts, vst.flavor='v2', batch_var='batch',
            latent_var_nonreg='percent', cell_attr=cell_attr,
            n_cells=ncol(counts), n_genes=nrow(counts), min_cells=5,
            return_gene_attr=TRUE, return_corrected_umi=TRUE, verbosity=0),
  stop('unsupported design-reference mode')
)
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
write_f64 <- function(name, value, matrix=FALSE) {
  con <- file(file.path(args[2], paste0(name, '.f64le')), 'wb')
  writeBin(if (matrix) as.double(t(value)) else as.double(value),
           con, size=8L, endian='little')
  close(con)
}
writeLines(rownames(result$model_pars), file.path(args[2], 'fit_genes.txt'))
writeLines(rownames(result$model_pars_fit), file.path(args[2], 'genes.txt'))
writeLines(colnames(result$model_pars), file.path(args[2], 'columns.txt'))
write_f64('model_pars', result$model_pars, TRUE)
write_f64('model_pars_fit', result$model_pars_fit, TRUE)
write_f64('residual', result$y, TRUE)
write_f64('residual_mean', result$gene_attr$residual_mean)
write_f64('residual_variance', result$gene_attr$residual_variance)
write_f64('outliers', result$model_pars_outliers)
if (length(result$model_pars_nonreg)) {
  writeLines(colnames(result$model_pars_nonreg),
             file.path(args[2], 'nonreg_columns.txt'))
  write_f64('model_pars_nonreg', result$model_pars_nonreg, TRUE)
}
if (!is.null(result$umi_corrected)) {
  con <- file(file.path(args[2], 'corrected_umi.i32le'), 'wb')
  writeBin(as.integer(t(as.matrix(result$umi_corrected))), con,
           size=4L, endian='little')
  close(con)
}
