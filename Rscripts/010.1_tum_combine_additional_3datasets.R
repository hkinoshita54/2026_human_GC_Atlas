# combine tumor cell objects from various data sets for NMF
# 2026-05-15

# load packages ----
library(tidyverse)
library(Seurat)
library(Matrix)

# load data ----
seu_list <- list()

## Tsutsumi et al. GSE201347
seu <- readRDS("../2026_Tsutsumi_GSE201347/RDSfiles/seu_020_tum.RDS")
str(seu[[]])
seu@meta.data <- seu@meta.data %>% select(
  orig.ident,
  nCount_RNA,
  nFeature_RNA,
  sample = gsm,
  patient,
  Lauren,
  percent.mt,
  S.Score,
  G2M.Score,
  Phase,
  EBV
)
seu$dataset <- "GSE201347"
str(seu[[]])
seu_list[[1]] <- seu

## Gao et al. GSE270680
seu <- readRDS("../2026_Gao_GSE270680/RDSfiles/seu_000_tum.RDS")
str(seu[[]])
seu@meta.data <- seu@meta.data %>% select(
  orig.ident = library,
  nCount_RNA,
  nFeature_RNA,
  sample = library,
  patient,
  Lauren = Subtype,
  percent.mt = percent.mito,
  S.Score,
  G2M.Score,
  Phase,
  celltype2 = subCluster
)
seu@meta.data$Lauren <- case_when(
  seu@meta.data$Lauren == "Intestinal" ~ "Intestinal",
  seu@meta.data$Lauren == "Indeterminate" ~ "Mixed",
  seu@meta.data$Lauren == "Diffuse" ~ "Diffuse"
)
seu$dataset <- "GSE270680"
str(seu[[]])
seu_list[[2]] <- seu

## Zhao et al. OMIX013242
seu <- readRDS("../2026_Zhao_OMIX013242/RDSfiles/seu_020_epi.RDS")
str(seu[[]])
seu@meta.data <- seu@meta.data %>% select(
  orig.ident,
  nCount_RNA,
  nFeature_RNA,
  sample = sample_id,
  patient,
  ACRG,
  TCGA,
  percent.mt,
  S.Score,
  G2M.Score,
  Phase
)
seu$dataset <- "OMIX013242"
str(seu[[]])
seu_list[[3]] <- seu

## merge data sets
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
str(seu[[]])
seu <- JoinLayers(seu)
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$data.2 <- NULL
saveRDS(seu, "RDSfiles/seu_010.1_tum2.RDS")
