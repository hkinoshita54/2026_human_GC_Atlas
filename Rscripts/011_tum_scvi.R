# Integrate subset of data sets
# 2026-03-30

# Setting up ----

## Make directories
analysis_step <- "011_tum_scvi"
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## helper
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load packages ----
library(tidyverse)
library(Seurat)
source("Rscripts/helpers.R")

# load data ----
seu <- readRDS("RDSfiles/seu_010_tum.RDS")
latent <- read.csv("adata/tum/scvi_latent.csv", row.names = 1, check.names = FALSE)

# add scvi latent into seurat ----
seu[["scvi"]] <- CreateDimReducObject(
  embeddings = as.matrix(latent),
  key = "scVI_",
  assay = DefaultAssay(seu)
)

# cluster using scvi latent ----
seu <- NormalizeData(seu)

seu <- FindNeighbors(seu, reduction = "scvi", dims = 1:ncol(latent))
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, reduction = "scvi", dims = 1:ncol(latent), reduction.name = "umap.scvi")

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

add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "IGKC", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF","HBB")
add_feat <- c("IGKC")
sapply(add_feat, save_fp, seu, fp_path)

# seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

VlnPlot(seu, features = c("nFeature_RNA"), pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("vln_nFeature.pdf", path = plot_path, width = 150, height = 35, units = "mm")

# re-cluster after removing Prietal, Tuft, and EEC, and Plasma cell contamination ----
seu <- subset(seu, seurat_clusters %in% c(6,20,23,24), invert = TRUE)
