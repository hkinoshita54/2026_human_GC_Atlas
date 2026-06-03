# Continued from 201.1_...
# 2026-05-25
# characterize MPs by enrichment test (hypergeometric test)

# Settings ----
## Make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "207_metaprogram_PB"

plot_path <- file.path(wd, "plot", analysis_step)
res_path <- file.path(wd, "result", analysis_step)
fp_path <- file.path(wd, "plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
library(DESeq2)

# pseudobulk ----
seu <- readRDS("RDSfiles/seu_203_tum_combined.RDS")
pb_counts <- AggregateExpression(seu, group.by = "sample", assays = "RNA", slot = "counts")$RNA # be sure to use raw counts

## extract sample meta data from Seurat object
sample_meta <- seu@meta.data %>%
  group_by(sample) %>%
  summarise(
    dataset = dplyr::first(dataset),
    patient = dplyr::first(patient),
    Lauren = dplyr::first(Lauren),
    TCGA = dplyr::first(TCGA),
    ACRG = dplyr::first(ACRG),
    layer = dplyr::first(layer),
    diagnosis = dplyr::first(diagnosis),
    EBV = dplyr::first(EBV),
    n_cells = n()
  ) %>%
  column_to_rownames("sample")
sample_meta <- sample_meta[colnames(pb_counts), , drop = FALSE]

## normalize counts by using DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(pb_counts), colData   = sample_meta, design    = ~ 1)
vsd <- vst(dds, blind = TRUE)
pb_vst <- assay(vsd)

# compute MP scores as mean z-scores of the MP genes ----
mp_res <- readRDS("RDSfiles/meta_program_results_2.rds")
mp_gene_lists <- mp_res$MP_list

## helper
score_gene_set <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  colMeans(expr[genes, , drop = FALSE])
}

## compute MP scores
pb_z <- t(scale(t(pb_vst)))
mp_scores <- sapply(mp_gene_lists, function(g) {
  score_gene_set(pb_z, g)
})
mp_scores <- as.data.frame(mp_scores)

## add to sample meta data
all(rownames(sample_meta) == rownames(mp_scores))
dat_all <- bind_cols(sample_meta, mp_scores)

# Diffuse vs Intestinal ----
dat <- dat_all %>% 
  filter(Lauren %in% c("Diffuse", "Intestinal"))
dat$Lauren <- factor(dat$Lauren, levels = c("Intestinal", "Diffuse"))
results <- lapply(colnames(mp_scores), function(mp) {
  fit <- lm(
    as.formula(paste0(mp, " ~ Lauren + dataset")),
    data = dat
  )
  broom::tidy(fit) %>%
    filter(grepl("Lauren", term)) %>%
    mutate(MP = mp)
}) %>% bind_rows()
openxlsx2::write_xlsx(results, file.path(res_path, "Lauren_MP.xlsx"))

# Deep vs Superficial ----
dat <- dat_all %>% 
  filter(layer %in% c("D", "S"))
dat$layer <- factor(dat$layer, levels = c("S", "D"))
results <- lapply(colnames(mp_scores), function(mp) {
  fit <- lm(
    as.formula(paste0(mp, " ~ layer")),
    data = dat
  )
  broom::tidy(fit) %>%
    filter(grepl("layer", term)) %>%
    mutate(MP = mp)
}) %>% bind_rows()
openxlsx2::write_xlsx(results, file.path(res_path, "layer_MP.xlsx"))

# early vs advanced ----
dat <- dat_all %>% 
  filter(!is.na(diagnosis))
dat$diagnosis <- factor(dat$diagnosis, levels = c("AGC", "EGC"))
results <- lapply(colnames(mp_scores), function(mp) {
  fit <- lm(
    as.formula(paste0(mp, " ~ diagnosis")),
    data = dat
  )
  broom::tidy(fit) %>%
    filter(grepl("diagnosis", term)) %>%
    mutate(MP = mp)
}) %>% bind_rows()
openxlsx2::write_xlsx(results, file.path(res_path, "advanced_early_MP.xlsx"))

# TCGA ----
dat <- dat_all %>% 
  filter(!is.na(TCGA))
dat$TCGA <- factor(dat$TCGA, levels = c("CIN", "GS", "MSI", "EBV"))
results <- lapply(colnames(mp_scores), function(mp) {
  fit <- lm(
    as.formula(paste0(mp, " ~ TCGA + dataset")),
    data = dat
  )
  broom::tidy(fit) %>%
    filter(grepl("TCGA", term)) %>%
    mutate(MP = mp)
}) %>% bind_rows()
openxlsx2::write_xlsx(results, file.path(res_path, "TCGA_MP.xlsx"))

# ACRG ----
dat <- dat_all %>% 
  filter(!is.na(ACRG))
dat$ACRG <- factor(dat$ACRG, levels = c("MSS/TP53-", "MSS/TP53+", "EMT", "MSI"))
results <- lapply(colnames(mp_scores), function(mp) {
  fit <- lm(
    as.formula(paste0(mp, " ~ ACRG + dataset")),
    data = dat
  )
  broom::tidy(fit) %>%
    filter(grepl("ACRG", term)) %>%
    mutate(MP = mp)
}) %>% bind_rows()
openxlsx2::write_xlsx(results, file.path(res_path, "ACRG_MP.xlsx"))