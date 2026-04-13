# preparation for convert seurat object to anndata

# load packages ----
library(tidyverse)
library(Seurat)
library(Matrix)
source("Rscripts/harmonize_symbols.R")
source("Rscripts/helpers.R")

# output directory
out_dir <- "out_data/epi"
fs::dir_create(out_dir)

# load data ----
seu_list <- list()

## Kang et al. GSE206785
seu <- readRDS("../2025_Kang/RDSfiles/seu_030_epi.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$celltype1 <- NULL
seu$seurat_clusters <- NULL
str(seu[[]])
seu_list[[1]] <- seu

## Kumar et al. GSE183904
seu <- readRDS("../2025_Kumar_2/RDSfiles/seu_030_epi.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$celltype1 <- NULL
seu$seurat_clusters <- NULL
str(seu[[]])
seu_list[[2]] <- seu

## Jeong et al. GSE167297
seu <- readRDS("../2023_Jeong/RDSfiles/seu_020_epi.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$celltype1 <- NULL
seu$celltype2 <- NULL
seu$RNA_snn_res.1 <- NULL
seu$seurat_clusters <- NULL
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$Lauren <- "Diffuse"
str(seu[[]])
seu_list[[3]] <- seu

## Jiang et al. GSE163558
seu <- readRDS("../2026_Jiang_GSE163558/RDSfiles/seu_020_epi.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$sample <- seu$orig.ident
seu$celltype1 <- NULL
seu$seurat_clusters <- NULL
str(seu[[]])
seu_list[[4]] <- seu

## Kim et al. GSE150290
seu <- readRDS("../2026_Kim/RDSfiles/seu_020_epi.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$celltype1 <- NULL
seu$seurat_clusters <- NULL
seu$RNA_snn_res.1 <- NULL
str(seu[[]])
seu_list[[5]] <- seu

## Sun et al. OMIX001073
seu <- readRDS("../2026_Sun_OMIX001073/RDSfiles/seu_020_epi.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = batch)
seu$n_counts <- NULL
seu$TCR_BCR <- NULL
seu$cluster <- NULL
seu$cluster2 <- NULL
seu$celltype1 <- NULL
seu$seurat_clusters <- NULL
str(seu[[]])
seu_list[[6]] <- seu

## merge 6 data sets
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
str(seu[[]])
seu <- JoinLayers(seu)

## adjust meta data
levels(seu$tissue_type)
seu$tissue_type[seu$tissue_type == "NT"] <- "N"
seu$tissue_type <- droplevels(seu$tissue_type)
table(seu$tissue_type, seu$dataset)

saveRDS(seu, "RDSfiles/seu_020_epi.RDS")

# convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/

# save metadata table
seu$barcode <- colnames(seu)
# seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
# seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file = file.path(out_dir, "seu_metadata.csv"), quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
writeMM(counts_matrix, file = file.path(out_dir, 'seu_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
# write.csv(seu@reductions$pca@cell.embeddings, file = file.path(out_dir, 'seu_pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = file.path(out_dir, 'seu_gene_names.csv'),
  quote = F, row.names = F, col.names = F
)

