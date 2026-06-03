# Continued from 201.1_...
# 2026-05-22
# characterize MPs by enrichment test (hypergeometric test)

# Settings ----
## Make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "205_metaprogram_annotation"

plot_path <- file.path(wd, "plot", analysis_step)
res_path <- file.path(wd, "result", analysis_step)
fp_path <- file.path(wd, "plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
library(msigdbr)

# load data ----
mp_res <- readRDS("RDSfiles/meta_program_results_2.rds")
seu <- readRDS("RDSfiles/seu_203_tum_combined.RDS")

## extract mp genes
mp_gene_sets <- mp_res$MP_list

## define universal background
counts <- LayerData(seu, assay = "RNA", layer = "counts")
gene_detect_frac <- Matrix::rowMeans(counts > 0)
bg_genes <- names(gene_detect_frac)[gene_detect_frac >= 0.01]
length(bg_genes)

## MSigDB collections for enrichment test
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
msig_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msig_c8 <- msigdbr(species = "Homo sapiens", category = "C8")
msig_all <- bind_rows(msig_h, msig_gobp, msig_c8) %>%
  select(gs_cat, gs_subcat, gs_name, gene_symbol) %>%
  distinct() %>%
  filter(gene_symbol %in% bg_genes)

# enrichment test (hypergeometric test) ----
## helper
run_hypergeom_enrichment <- function(mp_genes, genesets, bg_genes) {
  mp_genes <- intersect(unique(mp_genes), bg_genes)
  N <- length(unique(bg_genes))
  n <- length(mp_genes)
  
  genesets %>%
    group_by(gs_cat, gs_subcat, gs_name) %>%
    summarise(
      geneset_genes = list(unique(gene_symbol)),
      geneset_size = length(unique(gene_symbol)),
      overlap_genes = list(intersect(mp_genes, unique(gene_symbol))),
      overlap_size = length(intersect(mp_genes, unique(gene_symbol))),
      .groups = "drop"
    ) %>%
    filter(overlap_size > 0) %>%
    mutate(
      pval = phyper(
        q = overlap_size - 1,
        m = geneset_size,
        n = N - geneset_size,
        k = n,
        lower.tail = FALSE
      ),
      FDR = p.adjust(pval, method = "BH"),
      overlap_genes = map_chr(overlap_genes, ~ paste(.x, collapse = ", "))
    ) %>%
    arrange(FDR, desc(overlap_size))
}

## run enrichment test
mp_enrichment <- imap_dfr(
  mp_gene_sets,
  ~ run_hypergeom_enrichment(.x, msig_all, bg_genes) %>%
    mutate(MP = .y, .before = 1)
)
mp_enrichment_sig <- mp_enrichment %>%
  filter(FDR < 0.05)
openxlsx2::write_xlsx(mp_enrichment_sig, file.path(res_path, "mp_enrichment_sig.xlsx"))


# MP to MP correltation ----
mp_cols <- grep("^MP_[0-9]+_UCell$", colnames(seu@meta.data), value = TRUE)

mp_cor <- cor(
  seu@meta.data[, mp_cols],
  method = "spearman",
  use = "pairwise.complete.obs"
)

library(ComplexHeatmap)
library(circlize)
library(scico)

## rename
rownames(mp_cor) <- gsub("_UCell", "", rownames(mp_cor))
colnames(mp_cor) <- gsub("_UCell", "", colnames(mp_cor))

## color function
seq_cols <- scico(9, palette = "vik")
col_fun <- colorRamp2(seq(-1, 1, length.out = 9), seq_cols)

## heatmap
ht <- Heatmap(
  mp_cor,
  name = "Spearman\nrho",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title_position = "topcenter",
    legend_direction = "horizontal"
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", mp_cor[i, j]),
      x, y,
      gp = gpar(fontsize = 5)
    )
  }
)

pdf(file.path(plot_path, "MP_correlation_heatmap.pdf"), width = 4, height = 4)
draw(ht, heatmap_legend_side = "top", padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# MP enrichment by lineage/NOS clusters
library(presto)
mp_markers <- wilcoxauc(
  t(as.matrix(seu@meta.data[, mp_cols])),
  seu@meta.data$celltype1
)
write_csv(mp_markers, file = file.path(res_path, "mp_markers_annotated.csv"))
openxlsx2::write_xlsx(mp_markers, file.path(res_path, "mp_markers_annotated.xlsx"))
