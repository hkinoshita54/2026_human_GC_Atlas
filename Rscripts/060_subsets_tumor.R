# Integrate subset of data sets
# 2026-04-09

# Setting up ----

## Make directories
analysis_step <- "060_subsets_tumor"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(openxlsx2)
library(Seurat)
source("Rscripts/helpers.R")

## helper
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 4) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
seu_all <- readRDS("RDSfiles/seu_050.1_tum_harmony_patient.RDS")
seu <- subset(seu_all, subset = celltype2 == "Gastric")

## filter out patient with less than 100 cells
# patient_counts <- table(seu$patient)
# patients_keep <- names(patient_counts[patient_counts >= 100])
# seu <- subset(seu, subset = patient %in% patients_keep)
# seu$patient <- droplevels(seu$patient)
# table(seu$patient)

# Clustering w/ harmony per data set (1st) ----
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
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
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("cluster.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "copykat", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "copykat")
ggsave("copykat.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "dataset", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "dataset")
ggsave("dataset.pdf", path = plot_path, width = 40, height = 30, units = "mm")

add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("TME_score")
sapply(add_feat, save_fp, seu, fp_path)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))