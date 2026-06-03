# Continued from 208_... 
# 2026-05-30
# Create Seurat object for cellchat

# Settings ----
## Make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "311.1_cellchat_MP5_resid"

plot_path <- file.path(wd, "plot", analysis_step)
res_path <- file.path(wd, "result", analysis_step)
fp_path <- file.path(wd, "plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
library(CellChat)
source(file.path(ws, "common_Rscripts/helpers.R"))

# load seurat object ----
tum <- readRDS("RDSfiles/seu_203_tum_combined.RDS")
tum$ccgroup <- tum$MP5_resid_q

fib <- readRDS("RDSfiles/seu_111_fibro_tumor.RDS")
fib$ccgroup <- fib$celltype2

ec <- readRDS("RDSfiles/seu_121_ec_tumor.RDS")
ec$ccgroup <- ec$celltype2

myo <- readRDS("RDSfiles/seu_122_myocytes_tumor.RDS")
myo$ccgroup <- myo$celltype2

mye <- readRDS("RDSfiles/seu_131_mye_tumor.RDS")
mye$ccgroup <- mye$celltype2

seu <- merge(tum, list(fib, ec, myo, mye))
seu <- JoinLayers(seu)

## filter by ccgroup
ccgroup_levels <- c(unique(tum$MP5_resid_q), levels(fib$celltype2), levels(ec$celltype2), levels(myo$celltype2), levels(mye$celltype2))
seu$ccgroup <- factor(seu$ccgroup, levels = ccgroup_levels)
seu <- subset(seu, subset = ccgroup %in% c("LEC", "Stress_Mye", "Unclear_Mye"), invert = TRUE)

## keep sample_group with >=30 cells 
seu$n_sample_group <- seu@meta.data %>% 
  add_count(sample, ccgroup, name = "n_sample_group") %>% 
  pull(n_sample_group)
seu <- subset(seu, subset = n_sample_group >= 30)

## keep samaples with >= groups
samples_keep <- seu@meta.data %>%
  distinct(sample, ccgroup) %>%
  count(sample, name = "n_groups") %>%
  filter(n_groups >= 3) %>%
  pull(sample)
seu <- subset(seu, subset = sample %in% samples_keep)
seu_all <- seu # use seu_all as input for cellchat

# cellchat ----
load(file.path(ws, "aux_data/CellChatDB.new.RData"))
CellChatDB.use <- subsetDB(CellChatDB.new, search = c("Secreted Signaling"), key = "annotation")

seu_all <- NormalizeData(seu_all) # need normalized count for cellchat
samples <- unique(seu_all$sample)

cellchat_list <- list()
for (s in samples) {
  message("Running CellChat: ", s)
  seu <- subset(seu_all, subset = sample == s)
  seu$ccgroup <- droplevels(seu$ccgroup)
  seu$samples <- factor(seu$sample)
  cellchat <- createCellChat(object = seu, group.by = "ccgroup", assay = "RNA")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  if (nrow(cellchat@LR$LRsig) == 0) {
    message("Skipping ", s, ": no overexpressed LR pairs")
    next
  }
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 30)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  if (!is.null(cellchat@netP$prob) &&
      length(dim(cellchat@netP$prob)) == 3 &&
      dim(cellchat@netP$prob)[3] > 0) {
    
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
  } else {
    message("Skipping centrality for ", s, ": no valid pathway-level network")
  }
  
  cellchat_list[[s]] <- cellchat
}

saveRDS(cellchat_list, "RDSfiles/cellchat_list_311.1_MP5_resid.rds")

# collect information from CellChat results ----

cell_counts_all <- map_dfr(names(cellchat_list), function(s) {
  cellchat <- cellchat_list[[s]]
  
  tibble(
    sample = s,
    ccgroup = levels(cellchat@idents),
    n_cells = as.integer(table(cellchat@idents)[levels(cellchat@idents)])
  )
})

testable_pairs_all <- map_dfr(names(cellchat_list), function(s) {
  cellchat <- cellchat_list[[s]]
  groups <- levels(cellchat@idents)
  
  expand_grid(
    sample = s,
    source = groups,
    target = groups
  )
})

lr_all <- map_dfr(names(cellchat_list), function(s) {
  
  cellchat <- cellchat_list[[s]]
  
  prob_arr <- cellchat@net$prob
  pval_arr <- cellchat@net$pval
  
  lr_tbl <- as_tibble(as.data.frame.table(prob_arr, responseName = "prob")) %>%
    rename(
      source = Var1,
      target = Var2,
      interaction_name = Var3
    ) %>%
    mutate(
      sample = s,
      pval = as.vector(pval_arr),
      significant = prob > 0 & pval < 0.05
    ) %>%
    left_join(
      cellchat@LR$LRsig %>%
        as_tibble() %>%
        select(interaction_name, pathway_name, ligand, receptor),
      by = "interaction_name"
    ) %>%
    relocate(
      sample, source, target,
      ligand, receptor, pathway_name, interaction_name,
      prob, pval, significant
    )
  
  lr_tbl
})

openxlsx2::write_xlsx(cell_counts_all, file.path(res_path, "cell_counts_per_sample_ccgroup.xlsx"))
openxlsx2::write_xlsx(testable_pairs_all, file.path(res_path, "testable_source_target_pairs.xlsx"))
# openxlsx2::write_xlsx(lr_all, file.path(res_path, "LR_all_per_sample.xlsx"))

# summarize across samples ----
lr_summary <- lr_all %>%
  group_by(source, target, ligand, receptor, pathway_name, interaction_name) %>%
  summarise(
    n_testable_samples = n_distinct(sample),
    n_significant_samples = n_distinct(sample[significant]),
    fraction_significant = n_significant_samples / n_testable_samples,
    median_prob_among_testable = median(prob, na.rm = TRUE),
    median_prob_among_significant = if_else(
      n_significant_samples > 0,
      median(prob[significant], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_significant_samples),
    desc(fraction_significant),
    desc(median_prob_among_significant)
  )

openxlsx2::write_xlsx(lr_summary, file.path(res_path, "LR_summary_across_samples.xlsx"))

# MP5resid_high as target ----
mp5_target_compare <- lr_summary %>%
  filter(target %in% c("MP5_4", "MP5_1_2")) %>%
  select(
    source, target, ligand, receptor, pathway_name, interaction_name,
    n_testable_samples, n_significant_samples, fraction_significant,
    median_prob_among_testable, median_prob_among_significant
  ) %>%
  pivot_wider(
    names_from = target,
    values_from = c(
      n_testable_samples,
      n_significant_samples,
      fraction_significant,
      median_prob_among_testable,
      median_prob_among_significant
    )
  ) %>%
  mutate(
    delta_fraction_significant =
      fraction_significant_MP5_4 - fraction_significant_MP5_1_2,
    delta_median_prob_among_testable =
      median_prob_among_testable_MP5_4 - median_prob_among_testable_MP5_1_2,
    delta_median_prob_among_significant =
      median_prob_among_significant_MP5_4 - median_prob_among_significant_MP5_1_2
  ) %>%
  arrange(
    desc(delta_fraction_significant),
    desc(delta_median_prob_among_testable)
  )

openxlsx2::write_xlsx(mp5_target_compare, file.path(res_path, "MP5_4_vs_MP5_1_2_target_comparison.xlsx"))

# MP5resid_high as source ----
mp5_source_compare <- lr_summary %>%
  filter(source %in% c("MP5_4", "MP5_1_2")) %>%
  select(
    source, target, ligand, receptor, pathway_name, interaction_name,
    n_testable_samples, n_significant_samples, fraction_significant,
    median_prob_among_testable, median_prob_among_significant
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = c(
      n_testable_samples,
      n_significant_samples,
      fraction_significant,
      median_prob_among_testable,
      median_prob_among_significant
    )
  ) %>%
  mutate(
    delta_fraction_significant =
      fraction_significant_MP5_4 - fraction_significant_MP5_1_2,
    delta_median_prob_among_testable =
      median_prob_among_testable_MP5_4 - median_prob_among_testable_MP5_1_2,
    delta_median_prob_among_significant =
      median_prob_among_significant_MP5_4 - median_prob_among_significant_MP5_1_2
  ) %>%
  arrange(
    desc(delta_fraction_significant),
    desc(delta_median_prob_among_testable)
  )

openxlsx2::write_xlsx(mp5_source_compare, file.path(res_path, "MP5_4_vs_MP5_1_2_source_comparison.xlsx"))

# correlation with sequencing depth ----
library(ggpubr)
Idents(seu_all) <- "ccgroup"
FeatureScatter(seu_all, "MP5_resid", "nFeature_RNA", pt.size = 1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_MP5resid_nFeature.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

VlnPlot(seu_all, features = "nFeature_RNA", group.by = "MP5_resid_state", pt.size = 0)  &
  theme_panel() & NoLegend() &
  labs(title = "nFeature_RNA", x = NULL, y = "nFeature_RNA") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("vln_MP5resid_nFeature.pdf", path = plot_path, width = 30, height = 35, units = "mm", device = cairo_pdf)
