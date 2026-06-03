# 2026-05-27

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "122_myocytes_tumor"

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
seu <- subset(seu, subset = celltype1 %in% c("Myo") & tissue_type == "T")

# Clustering w/ harmony per patient (1st) ----
## remove patients with less than 10 cells
df <- table(seu$patient) %>% as.data.frame()
keep <- df$Var1[df$Freq >= 10]
seu <- subset(seu, subset = patient %in% keep)
seu$sample <- droplevels(seu$sample)

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
seu <- FindClusters(seu, resolution = 0.02, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("cluster_res0.02.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

# features <- readLines(file.path(ws, "aux_data/gene_set/global_markers.txt"))
# sapply(features, save_fp, seu, fp_path, pt.size = 4)

add_feat <- c("RGS5", "CSPG4", "MCAM", "PDGFRB", "DES", "ACTA2", "TAGLN", "MYH11", "CNN1", "TPM2", "COL1A1", "COL3A1", "ANGPT2", "CXCL12")
sapply(add_feat, save_fp, seu, fp_path, pt.size = 4)

add_feat <- c("PECAM1", "PTPRC", "CD3E", "CD79", "TYROBP", "PDGFRA", "POSTN")
sapply(add_feat, save_fp, seu, fp_path, pt.size = 4)

add_feat <- c("nFeature_RNA")
sapply(add_feat, save_fp, seu, fp_path, pt.size = 4)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_res0.2.png", path = plot_path, width = 6, height = 3, units = "in", dpi = 150)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.02.csv"))

# 2nd clustering ----
seu1 <- seu
seu <- subset(seu1, subset = seurat_clusters %in% c(3, 4), invert = TRUE) # remove low fibro and ec contamination

seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
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
seu <- FindClusters(seu, resolution = 0.02, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("cluster_2_res0.02.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt_2.pdf", path = plot_path, width = 25, height = 30, units = "mm")

fp_path_2 <- file.path(plot_path, "feature_plot_2")
fs::dir_create(fp_path_2)
# features <- readLines(file.path(ws, "aux_data/gene_set/global_markers.txt"))
# sapply(features, save_fp, seu, fp_path_2, pt.size = 4)

add_feat <- c("RGS5", "CSPG4", "MCAM", "PDGFRB", "DES", "ACTA2", "TAGLN", "MYH11", "CNN1", "TPM2", "COL1A1", "COL3A1", "ANGPT2", "CXCL12")
sapply(add_feat, save_fp, seu, fp_path_2, pt.size = 4)

add_feat <- c("PECAM1", "PTPRC", "CD3E", "CD79", "TYROBP", "PDGFRA", "POSTN")
sapply(add_feat, save_fp, seu, fp_path_2, pt.size = 4)

add_feat <- c("nFeature_RNA")
sapply(add_feat, save_fp, seu, fp_path_2, pt.size = 4)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_2_res0.01.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_2_res0.02.csv"))

# Add celltype2 annotation ----
seu2 <- seu
# seu <- subset(seu, subset = seurat_clusters %in% c(), invert = TRUE)
seu$celltype2 <- NA_character_
seu$celltype2[seu$seurat_clusters %in% c(0,3)] <- "vSMC" 
seu$celltype2[seu$seurat_clusters %in% c(1,2)] <- "Peri" 
seu$celltype2 <- factor(seu$celltype2, levels = c("vSMC", "Peri"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))

features <- c(
  "MYH11",
  "ACTA2",
  "DCN", 
  "RGS5"
)
DotPlot(seu, group.by = "celltype2", features = features, dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(title = "% Expr", title.position = "top"),
    colour = guide_colorbar(title = "Avg Expr", title.position = "top", barheight = grid::unit(10, "mm"), barwidth  = grid::unit(3, "mm"))
  ) +
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 40, height = 33, units = "mm")

# save ----
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
keep_cols <- grep("RNA_snn_res.", names(seu[[]]), invert = TRUE)
seu@meta.data <- seu@meta.data[, keep_cols]
saveRDS(seu, file = "RDSfiles/seu_122_myocytes_tumor.RDS")

add_feat <- c("ACTG2", "MYLK")
sapply(add_feat, save_fp, seu2, fp_path_2, pt.size = 4)