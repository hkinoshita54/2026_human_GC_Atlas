# 2026-05-28

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "131_mye_tumor"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source(file.path(ws, "common_Rscripts/helpers.R"))
source(file.path(ws, "common_Rscripts/harmonize_symbols.R"))

# Load data ----
## combine all the data sets
seu_list <- list()

## Kang et al. GSE206785
seu <- readRDS("../2025_Kang/RDSfiles/seu_050_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Kumar et al. GSE183904
seu <- readRDS("../2025_Kumar_2/RDSfiles/seu_050_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Jeong et al. GSE167297
seu <- readRDS("../2023_Jeong/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
seu$Lauren <- "Diffuse"
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Jiang et al. GSE163558
seu <- readRDS("../2026_Jiang_GSE163558/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Kim et al. GSE150290
seu <- readRDS("../2026_Kim/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Sun et al. OMIX001073
seu <- readRDS("../2026_Sun_OMIX001073/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Zhao et al. OMIX013242
seu <- readRDS("../2026_Zhao_OMIX013242/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Gao et al. GSE270680
seu <- readRDS("../2026_Gao_GSE270680/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Tsutsumi et al. GSE201347
seu <- readRDS("../2026_Tsutsumi_GSE201347/RDSfiles/seu_040_mye.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## merge 9 data sets
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
str(seu[[]])
seu$TCGA[seu$TCGA == "NA"] <- NA
seu$TCGA <- droplevels(seu$TCGA)
seu$ACRG <- case_when(
  seu$ACRG %in% c("MSS/TP53+", "MSS/TP53pos") ~ "MSS/TP53+",
  seu$ACRG %in% c("MSS/TP53-", "MSS/TP53neg") ~ "MSS/TP53-",
  seu$ACRG %in% c("EMT", "MSS/EMT") ~ "EMT",
  seu$ACRG %in% c("MSI") ~ "MSI",
)
seu$tissue_type <- case_when(
  seu$tissue_type %in% c("T") ~ "T",
  seu$tissue_type %in% c("N", "NT") ~ "N",
  TRUE ~ "T",
)
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
seu <- JoinLayers(seu)

# Clustering w/ harmony per patient (1st) ----
## only samples from tumor tissues
seu <- subset(seu, subset = tissue_type == "T")

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
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster_res0.2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

# features <- readLines(file.path(ws, "aux_data/gene_set/global_markers.txt"))
# sapply(features, save_fp, seu, fp_path)

add_feat <- c("TYROBP", "FCER1G", "C1QC", "APOE", "TREM2", "INHBA", "IL1B", "CD1C", "CLEC10A", "LAMP3", "CCR7", "MMP9", "S100A9", "CSF3R", "CSF1R", , "FSCN1", "LAMP3")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("CD3D", "CD3E", "KRT5", "KRT8", "KRT19", "EPCAM", "PECAM1", "VWF", "COL1A1", "DCN")
sapply(add_feat, save_fp, seu, fp_path)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_res0.2.png", path = plot_path, width = 10, height = 3, units = "in", dpi = 150)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.2.csv"))

# 2nd clustering ----
seu1 <- seu
seu <- subset(seu1, subset = seurat_clusters %in% c(3,12), invert = TRUE) # remove low quality cluster (low nFeature_RNA)

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
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster_2_res0.2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt_2.pdf", path = plot_path, width = 25, height = 30, units = "mm")

fp_path_2 <- file.path(plot_path, "feature_plot_2")
fs::dir_create(fp_path_2)

add_feat <- c("TYROBP", "FCER1G", "C1QC", "APOE", "TREM2", "INHBA", "IL1B", "CD1C", "CLEC10A", "LAMP3", "CCR7", "MMP9", "S100A9", "CSF3R", "CSF1R", "FSCN1", "LAMP3")
sapply(add_feat, save_fp, seu, fp_path_2)

add_feat <- c("CD3D", "CD3E", "KRT5", "KRT8", "KRT19", "EPCAM", "PECAM1", "VWF", "COL1A1", "DCN")
sapply(add_feat, save_fp, seu, fp_path_2)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_2_res0.1.png", path = plot_path, width = 10, height = 3, units = "in", dpi = 150)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_2_res0.2.csv"))

# Add celltype2 annotation ----
seu2 <- seu
seu$celltype2 <- NA_character_
seu$celltype2[seu$seurat_clusters %in% c(0)] <- "TAN"
seu$celltype2[seu$seurat_clusters %in% c(1)] <- "TAM_APOE"
seu$celltype2[seu$seurat_clusters %in% c(2,13)] <- "Mono_FCN1"
seu$celltype2[seu$seurat_clusters %in% c(3,12)] <- "cDC2"
seu$celltype2[seu$seurat_clusters %in% c(4,5,10)] <- "Infl_TAM"
seu$celltype2[seu$seurat_clusters %in% c(6)] <- "DC_LAMP3"
seu$celltype2[seu$seurat_clusters %in% c(7)] <- "cDC1"
seu$celltype2[seu$seurat_clusters %in% c(8)] <- "Stress_Mye"
seu$celltype2[seu$seurat_clusters %in% c(9)] <- "pDC"
seu$celltype2[seu$seurat_clusters %in% c(11)] <- "Unclear_Mye"
seu$celltype2 <- factor(seu$celltype2, levels = c("Mono_FCN1", "Infl_TAM", "TAM_APOE", "cDC1", "cDC2", "DC_LAMP3", "pDC", "TAN", "Stress_Mye", "Unclear_Mye"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype2.pdf", path = plot_path, width = 60, height = 30, units = "mm")

# save ----
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
keep_cols <- grep("RNA_snn_res.", names(seu[[]]), invert = TRUE)
seu@meta.data <- seu@meta.data[, keep_cols]
saveRDS(seu, file = "RDSfiles/seu_131_mye_tumor.RDS")
