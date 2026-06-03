# 2026-05-18
# NMF based metaprogram discovery as in Gavish et al. nature 2023, more faithful version
# https://github.com/tiroshlab/3ca/tree/main/ITH_hallmarks
# 3 additional data sets from 010.1_

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "200.1_NMF"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
library(NMF)
library(doParallel)
library(Matrix)
source(file.path(ws, "common_Rscripts/helpers.R"))

# load data ----
## this seurat object is created in 010_tum_seurat_to_adata.R, gene symbols are harmonized already
seu <- readRDS("RDSfiles/seu_010.1_tum2.RDS")
seu_list <- SplitObject(seu, split.by = "sample")
sapply(seu_list, ncol) # check cell number of each sample
seu_list <- seu_list[sapply(seu_list, ncol) >= 10]
sapply(seu_list, ncol) # check cell number of each sample

# run nmf for all the patients ----
## function ----
run_nmf <- function(seu, patient, nfeatures = 7000,
                    ranks = 4:9, nrun = 10,
                    nmf_method = "snmf/r",
                    options = paste0("vP", ncores),
                    count_layer = "counts") {
  
  message("Running NMF for ", patient, " (", ncol(seu), " cells)")
  
  counts <- LayerData(seu, layer = count_layer)
  
  ## Remove genes with zero counts across this sample
  genes_detected <- Matrix::rowSums(counts) > 0
  counts <- counts[genes_detected, , drop = FALSE]
  
  ## CPM normalization
  lib_size <- Matrix::colSums(counts)
  
  if (any(lib_size == 0)) {
    stop("Some cells have zero library size in patient: ", patient)
  }
  
  cpm <- t(t(counts) / lib_size * 1e6)
  
  log_cpm <- log2(cpm / 10 + 1)
  
  ## Select top 7000 genes by mean expression before centering
  gene_mean <- Matrix::rowMeans(log_cpm)
  
  nfeatures_use <- min(nfeatures, length(gene_mean))
  genes_use <- names(sort(gene_mean, decreasing = TRUE))[seq_len(nfeatures_use)]
  
  mtx <- log_cpm[genes_use, , drop = FALSE]
  
  ## Center each gene across cells, Set negative values to zero
  gene_center <- Matrix::rowMeans(mtx)
  mtx <- mtx - gene_center
  mtx[mtx < 0] <- 0
  mtx <- as.matrix(mtx)
  
  ## safety checks
  stopifnot(all(mtx >= 0))
  stopifnot(!any(is.na(mtx)))
  stopifnot(!any(is.infinite(mtx)))
  
  ## remove genes that became all-zero after centering
  keep_genes <- rowSums(mtx) > 0
  mtx <- mtx[keep_genes, , drop = FALSE]
  
  message("Using ", nrow(mtx), " genes after centering/non-negative transform")
  
  start <- Sys.time()
  
  nmf_res <- nmf(
    mtx,
    rank = ranks,
    method = nmf_method,
    nrun = nrun,
    .options = options,
    .pbackend = NULL
  )
  
  end <- Sys.time()
  
  message(
    "Finished ", patient, " in ",
    round(difftime(end, start, units = "mins"), 2), " min"
  )
  
  W_list <- lapply(names(nmf_res$fit), function(k) {
    W <- basis(nmf_res$fit[[k]])
    colnames(W) <- paste(patient, k, seq_len(ncol(W)), sep = "_")
    W
  })
  
  W_all <- do.call(cbind, W_list)
  
  return(W_all)
}


## run all patients ----
ncores <- 10

cl <- makeCluster(ncores)
registerDoParallel(cl)
nmf.options(pbackend = NULL)

Genes_nmf_w_basis <- list()

for (patient in names(seu_list)) {
  Genes_nmf_w_basis[[patient]] <- run_nmf(
    seu = seu_list[[patient]],
    patient = patient
  )
}

stopCluster(cl)
registerDoSEQ()

saveRDS(
  Genes_nmf_w_basis,
  "RDSfiles/Genes_nmf_w_basis_log2CPM10_centered_rank4_9_nrun10_min10cells_additional.RDS"
)
