# Integrate subset of data sets
# 2026-04-03

# Setting up ----

## Make directories
analysis_step <- "050_epi_harmony_patient"
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
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# load data ----
seu <- readRDS("RDSfiles/seu_020_epi.RDS")

## filter out patient with less than 100 cells
patient_counts <- table(seu$patient)
patients_keep <- names(patient_counts[patient_counts >= 100])
seu <- subset(seu, subset = patient %in% patients_keep)
seu$patient <- droplevels(seu$patient)
table(seu$patient)

## add copykat results
tum <- readRDS("RDSfiles/seu_010_tum.RDS")
ck <- tum[[]] %>% select(copykat)
seu <- AddMetaData(seu, metadata = ck)
table(seu$copykat, useNA = "ifany")
str(seu[[]])

# Clustering w/ harmony per data set (1st) ----
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
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

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "tissue_type")
ggsave("tissue_type.pdf", path = plot_path, width = 30, height = 30, units = "mm")

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

# save seurat object for later visualization ----
seu_1 <- seu
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
# seu <- DietSeurat(seu, dimreducs = "umap")
saveRDS(seu, "RDSfiles/seu_050_epi_harmony_patient.RDS")

# remove clearly differentiated, non-malignant clusters and subset to T > define as tumor cells ----
## parietal: c11, endocrine: c8 & 12, tuft: c15, goblet/IM: c17, chief: c9, immune contamination: c14
seu <- subset(seu, subset = seurat_clusters %in% c(11, 8,12, 15, 17, 9, 14), invert = TRUE) 
seu <- subset(seu, subset = tissue_type == "T")

## filter out patient with less than 50 cells
patient_counts <- table(seu$patient)
patients_keep <- names(patient_counts[patient_counts >= 50])
seu <- subset(seu, subset = patient %in% patients_keep)
seu$patient <- droplevels(seu$patient)
table(seu$patient) # > 71 patients remain

# Clustering "tumor cells" w/ harmony per data set (2nd) ----
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)
# seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- ScaleData(seu, features = hvg) # without regressing out cell cycle scores

npcs <- 50
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
ggsave("cluster0.2.pdf", path = plot_path, width = 45, height = 30, units = "mm")

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
ggsave("diagnosis.pdf", path = plot_path, width = 30, height = 30, units = "mm")
add_feat <- c("MKI67", "TOP2A", "MUC5AC", "TFF1", "MUC6", "PGA4", "ATP4B", "CHGA", "TRPM5", "TFF3", "CEACAM5","CEACAM6", "MSLN")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("PTPRC", "CD3E", "CD79A", "JCHAIN", "TYROBP", "S100A9", "MS4A2", "COL1A1", "PDGFRA", "MYH11", "VWF")
add_feat <- c("FABP1", "FABP2", "CDX1", "CDX2", "CDH1")
sapply(add_feat, save_fp, seu, fp_path)

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_tum.csv"))

VlnPlot(seu, features = c("nFeature_RNA"), pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("vln_nFeature.pdf", path = plot_path, width = 150, height = 35, units = "mm")

# UCell scoring ----
library(msigdbr)
library(UCell)
library(ggpubr)
H <- msigdbr(species    = "Homo sapiens", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

yap_list <- readRDS("RDSfiles/yap_list_human.RDS")
yap_list <- yap_list[4:5]

seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

seu <- AddModuleScore_UCell(
  seu,
  features = yap_list,
  chunk.size = 1000,
  ncores = 8,
)

## make ucell scores into z-scores
### change names for brevity
colnames(seu@meta.data)[colnames(seu@meta.data) == "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell"] <- "EMT_UCell"
colnames(seu@meta.data)[colnames(seu@meta.data) == "INFLAMMATORY_RESPONSE_UCell"] <- "INFL_RESP_UCell"

ucell_cols <- grep("_UCell$", colnames(seu[[]]), value = TRUE)

for (cc in ucell_cols) {
  seu[[cc]] <- as.numeric(scale(seu@meta.data[[cc]]))
}

## create a composite score
seu[["TME_score"]] <- rowMeans(
  seu@meta.data[, c(
    "EMT_UCell",
    "HYPOXIA_UCell",
    "INFL_RESP_UCell"
  )],
  na.rm = TRUE
)

add_feat <- c("nFeature_RNA")
sapply(add_feat, save_fp, seu, fp_path)

# annotation ----
seu$celltype2 <- NA_character_
seu$celltype2[seu$seurat_clusters %in% c(0)] <- "Intestinal"
seu$celltype2[seu$seurat_clusters %in% c(1,6,10,7)] <- "Gastric"
seu$celltype2[seu$seurat_clusters %in% c(2,5,11)] <- "Ambiguous"
seu$celltype2[seu$seurat_clusters %in% c(3)] <- "Neck_like"
seu$celltype2[seu$seurat_clusters %in% c(4)] <- "Goblet_like"
seu$celltype2[seu$seurat_clusters %in% c(8,9,12,13)] <- "PtSpecific"
seu$celltype2 <- factor(seu$celltype2, levels = c("Intestinal", "Gastric", "Ambiguous", "Neck_like", "Goblet_like", "PtSpecific"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes()
ggsave("celltype2.pdf", path = plot_path, width = 40, height = 30, units = "mm")


# save tumor subset ----
seu_2 <- seu
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
# seu <- DietSeurat(seu, dimreducs = "umap")
saveRDS(seu, "RDSfiles/seu_050.1_tum_harmony_patient.RDS")

# convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/
library(Matrix)

# output directory
out_dir <- "out_data/tum_harmony_patient"
fs::dir_create(out_dir)

# save metadata table
seu$barcode <- colnames(seu)
# seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
# seu$UMAP_2 <- seu@reductions$umap@cell.embeddings[,2]
write.csv(seu@meta.data, file = file.path(out_dir, "seu_metadata.csv"), quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- LayerData(seu, assay = 'RNA', layer = 'counts')
Matrix::writeMM(counts_matrix, file = file.path(out_dir, 'seu_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
# write.csv(seu@reductions$pca@cell.embeddings, file = file.path(out_dir, 'seu_pca.csv'), quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene' = rownames(counts_matrix)), file = file.path(out_dir, 'seu_gene_names.csv'),
  quote = F, row.names = F, col.names = F
)
