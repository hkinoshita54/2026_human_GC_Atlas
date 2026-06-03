# 2026-05-26

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "110_fibro"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source(file.path(ws, "common_Rscripts/helpers.R"))

# Load data ----
seu <- readRDS("RDSfiles/seu_100_str.RDS")
seu <- subset(seu, subset = celltype1 == "Fibro")

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

# Clustering w/ harmony per patient (1st) ----
## remove patients with less than 10 cells
df <- table(seu$patient) %>% as.data.frame()
keep <- df$Var1[df$Freq >= 10]
seu <- subset(seu, subset = patient %in% keep)

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000)
hvg <- VariableFeatures(seu)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))

npcs <- 30
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("cluster_res0.2.pdf", path = plot_path, width = 50, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

features <- readLines(file.path(ws, "aux_data/gene_set/global_markers.txt"))
sapply(features, save_fp, seu, fp_path)

# add_feat <- c("")
# sapply(add_feat, save_fp, seu, fp_path)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.2.csv"))
# 
# 
# # Add celltype1 annotation ----
# seu1 <- seu
# seu <- subset(seu, subset = seurat_clusters %in% c(14), invert = TRUE) # ambient RNA with immune or doublets
# seu$celltype1 <- NA_character_
# seu$celltype1[seu$seurat_clusters %in% c(1,3,4,12,13,16,19,22,25)] <- "EC"
# seu$celltype1[seu$seurat_clusters %in% c(2,5,6,7,9,10,11,17,18,24)] <- "Fibro"
# seu$celltype1[seu$seurat_clusters %in% c(0,8,15,20,21)] <- "Myo"
# seu$celltype1[seu$seurat_clusters %in% c(23)] <- "Glia"
# seu$celltype1 <- factor(seu$celltype1, levels = c("EC", "Fibro", "Myo", "Glia"))
# DimPlot(seu, group.by = "celltype1", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
#   theme_panel() & NoAxes() &
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
# ggsave("celltype1.pdf", path = plot_path, width = 40, height = 30, units = "mm")
# 
# # save ----
# seu[["RNA"]]$data <- NULL
# seu[["RNA"]]$scale.data <- NULL
# seu$RNA_snn_res.0.2 <- NULL
# seu$RNA_snn_res.0.5 <- NULL
# # seu$RNA_snn_res.1 <- NULL
# saveRDS(seu, file = "RDSfiles/seu_100_str.RDS")
