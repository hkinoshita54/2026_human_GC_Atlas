# 2026-05-19

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "203_visualize_metaprogram"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source(file.path(ws, "common_Rscripts/helpers.R"))

# Load data ----
## combine all the data sets
seu <- readRDS("RDSfiles/seu_010_tum.RDS")
seu2 <- readRDS("RDSfiles/seu_010.1_tum2.RDS")
seu <- merge(x = seu, y = seu2)
rm(seu2)
seu <- JoinLayers(seu)

### organize meta.data
str(seu[[]])
seu@meta.data <- seu@meta.data %>% 
  select(
    orig.ident, nCount_RNA, nFeature_RNA, dataset, sample, patient, Lauren, TCGA, ACRG, 
    percent.mt, S.Score, G2M.Score, Phase, copykat, layer, diagnosis, EBV, celltype2
  )

seu$Lauren <- case_when(
  seu$Lauren %in% c("diffuse", "Diffuse") ~ "Diffuse",
  seu$Lauren %in% c("intestinal", "Intestinal") ~ "Intestinal",
  seu$Lauren %in% c("mixed", "Mixed") ~ "Mixed",
  seu$Lauren %in% c("Metastatic") ~ "Metastatic",
)

seu$TCGA[seu$TCGA == "NA"] <- NA_character_

seu$ACRG[seu$ACRG == "NA"] <- NA_character_
seu$ACRG <- case_when(
  seu$ACRG %in% c("MSS/TP53+", "MSS/TP53pos") ~ "MSS/TP53+",
  seu$ACRG %in% c("MSS/TP53-", "MSS/TP53neg") ~ "MSS/TP53-",
  seu$ACRG %in% c("EMT", "MSS/EMT") ~ "EMT",
  seu$ACRG %in% c("MSI") ~ "MSI",
)

seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)
str(seu[[]])

saveRDS(seu, "RDSfiles/seu_203_tum_combined.RDS")

# Clustering ----
## remove samples with less than 10 cells
df <- table(seu$sample) %>% as.data.frame()
keep <- df$Var1[df$Freq >= 10]
seu <- subset(seu, subset = sample %in% keep)

## harmony integration per patient
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
hvg <- VariableFeatures(seu)

seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))

npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Cluster") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

## Check markers by feature plots
features <- readLines(file.path(ws, "aux_data/gene_set/global_markers.txt"))
sapply(features, save_fp, seu, fp_path)

add_feat <- c("CHGA","CHGB","TRPM5","DCLK1","FABP1","FABP2","TFF3","MUC2")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("MKI67", "TOP2A")
sapply(add_feat, save_fp, seu, fp_path)

# Add celltype1 annotation ----
seu$celltype1 <- "NOS4"
seu$celltype1[seu$seurat_clusters %in% c(0,6)] <- "Pit"
seu$celltype1[seu$seurat_clusters %in% c(3,9)] <- "Neck"
seu$celltype1[seu$seurat_clusters %in% c(22)] <- "Pariet"
seu$celltype1[seu$seurat_clusters %in% c(8)] <- "EEC"
seu$celltype1[seu$seurat_clusters %in% c(7,15)] <- "Entero"
seu$celltype1[seu$seurat_clusters %in% c(18)] <- "Goblet"
seu$celltype1[seu$seurat_clusters %in% c(1,4,12,13,21)] <- "NOS1"
seu$celltype1[seu$seurat_clusters %in% c(2,14,16,19,20,23)] <- "NOS2"
seu$celltype1[seu$seurat_clusters %in% c(5,10)] <- "NOS3"
seu$celltype1 <- factor(seu$celltype1, levels = c("Pit", "Neck", "Pariet", "EEC", "Entero", "Goblet", 
                                                  "NOS1", "NOS2", "NOS3", "NOS4"))
DimPlot(seu, group.by = "celltype1", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("celltype1.pdf", path = plot_path, width = 45, height = 30, units = "mm")

Idents(seu) <- "celltype1"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))

# visualize MPs ----
MP_results <- readRDS("RDSfiles/meta_program_results_2.rds")
library(UCell)
library(msigdbr)
# library(ggpubr)

## UCell scoring
H <- msigdbr(species    = "Homo sapiens", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

yap_list <- readRDS("RDSfiles/yap_list_human.RDS")
yap_list <- yap_list[4:5]

MP_list <- MP_results$MP_list

seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

colnames(seu@meta.data)[colnames(seu@meta.data) == "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell"] <- "EMT_UCell"
colnames(seu@meta.data)[colnames(seu@meta.data) == "INFLAMMATORY_RESPONSE_UCell"] <- "INFL_RESP_UCell"

seu <- AddModuleScore_UCell(
  seu,
  features = yap_list,
  chunk.size = 1000,
  ncores = 8,
)

seu <- AddModuleScore_UCell(
  seu,
  features = MP_list,
  chunk.size = 1000,
  ncores = 8,
)

ucell_cols <- grep("_UCell$", colnames(seu[[]]), value = TRUE)

for (cc in ucell_cols) {
  seu[[cc]] <- as.numeric(scale(seu@meta.data[[cc]]))
}

sapply(ucell_cols, save_fp, seu, fp_path)

# save RDS ----
# seu <- DietSeurat(seu)
# seu[["RNA"]]$scale.data <- NULL
# seu[["RNA"]]$data <- NULL
saveRDS(seu, file = "RDSfiles/seu_203_tum_combined.RDS")

# convert Seurat object to anndata manually following the tutorial below ----
# https://smorabit.github.io/tutorials/8_velocyto/
out_dir <- "out_data/tum_9_datasets/"

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

# write gene names of MPs
mp_genes <- unique(unlist(MP_list))
mp_genes <- mp_genes[!is.na(mp_genes)]
write.table(
  data.frame('gene' = mp_genes), file = file.path(out_dir, 'MP_genes.csv'),
  quote = F, row.names = F, col.names = T
)
