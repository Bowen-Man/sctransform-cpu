args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('usage: Rscript r_runner.R DATA OUTPUT.json')
suppressPackageStartupMessages({library(Matrix); library(sctransform); library(jsonlite)})
data_dir <- args[1]
counts <- as(readMM(file.path(data_dir, 'counts.mtx')), 'CsparseMatrix')
rownames(counts) <- read.delim(file.path(data_dir, 'genes.tsv'))$gene
colnames(counts) <- read.delim(file.path(data_dir, 'cells.tsv'))$cell
set.seed(42L)
timing <- system.time(result <- vst(
  counts,
  vst.flavor='v2',
  n_cells=ncol(counts),
  n_genes=nrow(counts),
  min_cells=5,
  return_gene_attr=TRUE,
  return_corrected_umi=TRUE,
  verbosity=0
))
payload <- list(
  implementation='sctransform::vst',
  compute_wall_seconds=unname(timing[['elapsed']]),
  compute_cpu_seconds=unname(timing[['user.self']] + timing[['sys.self']]),
  retained_genes=nrow(result$model_pars_fit),
  fitted_genes=nrow(result$model_pars),
  fitted_cells=length(result$cells_step1),
  r_version=as.character(getRversion()),
  sctransform_version=as.character(packageVersion('sctransform')),
  glmgampoi_version=as.character(packageVersion('glmGamPoi'))
)
write_json(payload, args[2], auto_unbox=TRUE, pretty=TRUE)
