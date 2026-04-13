# Continued from 031_tum_sub_seurat_to_adata and pySCENIC output
# 2026-03-30

# Settings ----
## Make directories
wd <- getwd()
analysis_step <- "033_tum_sub_scenic"
plot_path <- file.path(wd, "plot", analysis_step)
res_path <- file.path(wd, "result", analysis_step)
fp_path <- file.path(wd, "plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source("Rscripts/helpers.R")
# library(ggrepel)

## helper
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
seu <- readRDS("RDSfiles/seu_030_tum_sub.RDS")
auc_mtx <- read_tsv(file = "adata/tum_sub/auc_mtx.txt") %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
auc_mtx <- auc_mtx[,Cells(seu)]
seu[["reg"]] <- CreateAssayObject(data = auc_mtx)
DefaultAssay(seu) <- "reg"

# clustering based on auc_mtx (w/o pca)
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
npcs <- 30
seu <- RunPCA(seu, features = rownames(seu), npcs = npcs, reduction.name = "pca_reg", verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:npcs, reduction = "pca_reg", reduction.name = "umap_scenic", verbose = FALSE)
seu <- FindNeighbors(seu, graph.name = "reg_nn", reduction = "pca_reg", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, graph.name = "reg_nn", resolution = 0.5, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.pdf", path = plot_path, width = 45, height = 30, units = "mm")

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

harmony_annotation <- readRDS("RDSfiles/030_harmony_annotation.RDS")
seu <- AddMetaData(seu, metadata = harmony_annotation)

DimPlot(seu, cols = "polychrome", group.by = "celltype2", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "harmony_annotation")
ggsave("harmony_annotation.pdf", path = plot_path, width = 40, height = 30, units = "mm")

# check regulons as markers ----
# DefaultAssay(seu) <- "reg"
# 
# add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
# sapply(add_feat, save_fp, seu, fp_path)
# 
# add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF")
# sapply(add_feat, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_reg.csv"))

# check RNA markers of scenic identified clusters ----
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)

add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF")
sapply(add_feat, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

