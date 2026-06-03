# 2026-06-01

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "151_b_plasma_tumor"

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
seu <- readRDS("../2025_Kang/RDSfiles/seu_070_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Kumar et al. GSE183904
seu <- readRDS("../2025_Kumar_2/RDSfiles/seu_070_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Jeong et al. GSE167297
seu <- readRDS("../2023_Jeong/RDSfiles/seu_060_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
seu$Lauren <- "Diffuse"
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Jiang et al. GSE163558
seu <- readRDS("../2026_Jiang_GSE163558/RDSfiles/seu_060_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Kim et al. GSE150290
seu <- readRDS("../2026_Kim/RDSfiles/seu_060_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Sun et al. OMIX001073
seu <- readRDS("../2026_Sun_OMIX001073/RDSfiles/seu_060_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Zhao et al. OMIX013242
seu <- readRDS("../2026_Zhao_OMIX013242/RDSfiles/seu_060_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
seu$seurat_clusters <- NULL
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Gao et al. GSE270680
seu <- readRDS("../2026_Gao_GSE270680/RDSfiles/seu_060_b_plasma.RDS")
seu <- harmonize_symbols(seu)
str(seu[[]])
dataset <- seu$dataset %>% unique() %>% as.character()
seu_list[[dataset]] <- seu

## Tsutsumi et al. GSE201347
seu <- readRDS("../2026_Tsutsumi_GSE201347/RDSfiles/seu_060_b_plasma.RDS")
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

## remove samples with less than 10 cells
df <- table(seu$sample) %>% as.data.frame()
keep <- df$Var1[df$Freq >= 10]
seu <- subset(seu, subset = sample %in% keep)

## filter genes by expression
keep <- Matrix::rowSums(seu[["RNA"]]$counts > 0) >= 20
seu <- subset(seu, features = rownames(seu)[keep])

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 2000)
hvg <- VariableFeatures(seu)
# seu <- ScaleData(seu, features = hvg)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))

npcs <- 30
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
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster_res0.5.pdf", path = plot_path, width = 45, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_res0.5.png", path = plot_path, width = 15, height = 3, units = "in", dpi = 150)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, max.cells.per.ident = 2000)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_res0.5.csv"))

add_feat <- c("IGHD", "GPR183", "IRAG2", "JCHAIN", "MKI67")
sapply(add_feat, save_fp, seu, fp_path)

# Add celltype2 annotation ----
# seu <- subset(seu2, subset = seurat_clusters %in% c(), invert = TRUE) # remove low quality cluster
seu$celltype2 <- NA_character_
seu$celltype2[seu$seurat_clusters %in% c(1)] <- "Bn_IGHD"
seu$celltype2[seu$seurat_clusters %in% c(0,5,8,10)] <- "Bm_GPR183"
seu$celltype2[seu$seurat_clusters %in% c(12)] <- "GCB_IRAG2"
seu$celltype2[seu$seurat_clusters %in% c(2,3,4,6,7,9,11,13,15)] <- "Plasma"
seu$celltype2[seu$seurat_clusters %in% c(14)] <- "Cycling_B"
seu$celltype2 <- factor(
  seu$celltype2,
  levels = c(
    "Bn_IGHD",
    "Bm_GPR183",
    "GCB_IRAG2",
    "Plasma",
    "Cycling_B"
  )
)
DimPlot(seu, group.by = "celltype2", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("celltype2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DotPlot(seu, group.by = "celltype2", features = add_feat, dot.scale = 2.3) + 
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
ggsave("dotplot.pdf", path = plot_path, width = 45, height = 35, units = "mm")

# save ----
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
keep_cols <- grep("RNA_snn_res.", names(seu[[]]), invert = TRUE)
seu@meta.data <- seu@meta.data[, keep_cols]
saveRDS(seu, file = "RDSfiles/seu_151_b_plasma_tumor.RDS")

seu1 <- seu[, !is.na(seu$annotation2)]
DimPlot(seu1, group.by = "annotation2", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("annotation2.pdf", path = plot_path, width = 40, height = 30, units = "mm")
