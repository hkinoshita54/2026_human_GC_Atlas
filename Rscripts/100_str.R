# Integrate subset of data sets
# 2026-04-07

# Setting up ----

## Make directories
analysis_step <- "100_str"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(openxlsx2)
library(Seurat)
source("Rscripts/harmonize_symbols.R")
source("Rscripts/helpers.R")

## helper
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
seu_list <- list()

## Kang et al. GSE206785
seu <- readRDS("../2025_Kang/RDSfiles/seu_040_str.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$seurat_clusters <- NULL
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Kumar et al. GSE183904
seu <- readRDS("../2025_Kumar_2/RDSfiles/seu_040_str.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Jeong et al. GSE167297
seu <- readRDS("../2023_Jeong/RDSfiles/seu_030_str.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$celltype2 <- NULL
seu$RNA_snn_res.1 <- NULL
seu$seurat_clusters <- NULL
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$Lauren <- "Diffuse"
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Jiang et al. GSE163558
seu <- readRDS("../2026_Jiang_GSE163558/RDSfiles/seu_030_str.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$sample <- seu$orig.ident
seu$seurat_clusters <- NULL
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Kim et al. GSE150290
seu <- readRDS("../2026_Kim/RDSfiles/seu_030_str.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = gsm)
seu$seurat_clusters <- NULL
seu$RNA_snn_res.1 <- NULL
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Sun et al. OMIX001073
seu <- readRDS("../2026_Sun_OMIX001073/RDSfiles/seu_030_str.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu@meta.data <- seu@meta.data %>% rename(sample = batch)
seu$n_counts <- NULL
seu$TCR_BCR <- NULL
seu$cluster <- NULL
seu$cluster2 <- NULL
seu$seurat_clusters <- NULL
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## merge 6 data sets
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
str(seu[[]])
seu <- JoinLayers(seu)

# Clustering w/o integration ----
# seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
# hvg <- VariableFeatures(seu)
# cc.genes <- Seurat::cc.genes
# seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
# seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
# npcs <- 50
# seu <- RunPCA(seu, npcs = npcs)
# seu <- FindNeighbors(seu, dims = 1:npcs)
# seu <- FindClusters(seu, resolution = 0.5)
# seu <- RunUMAP(seu, dims = 1:npcs)

# Clustering w/ harmony per data set (1st) ----
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$dataset)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))

npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("cluster.pdf", path = plot_path, width = 50, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "copykat", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "copykat")
ggsave("copykat.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "dataset", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "dataset")
ggsave("dataset.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "diagnosis", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "diagnosis")
ggsave("diagnosis.pdf", path = plot_path, width = 40, height = 30, units = "mm")

add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF")
sapply(add_feat, save_fp, seu, fp_path)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))


# re-clustering after removing obvious non-tumor cells (2nd) ----
seu <- subset(seu, subset = seurat_clusters %in% c(10,17,19,23,24,25), invert = TRUE)

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$dataset)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)
# seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- ScaleData(seu, features = hvg)

npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("cluster.pdf", path = plot_path, width = 50, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "copykat", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "copykat")
ggsave("copykat.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "dataset", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "dataset")
ggsave("dataset.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "diagnosis", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "diagnosis")
ggsave("diagnosis.pdf", path = plot_path, width = 40, height = 30, units = "mm")

add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF")
sapply(add_feat, save_fp, seu, fp_path)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers2.csv"))

# re-doing clustering without regressing out cc (3rd) ----
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$dataset)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)
# seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- ScaleData(seu, features = hvg)

npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("cluster.pdf", path = plot_path, width = 50, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "copykat", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "copykat")
ggsave("copykat.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "dataset", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "dataset")
ggsave("dataset.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "diagnosis", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "diagnosis")
ggsave("diagnosis.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "Platform", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Platform")
ggsave("Platform.pdf", path = plot_path, width = 40, height = 30, units = "mm")

add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF")
sapply(add_feat, save_fp, seu, fp_path)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers3.csv"))