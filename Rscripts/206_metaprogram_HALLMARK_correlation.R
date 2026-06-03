# Continued from 201.1_...
# 2026-05-24
# sample levels correlation between MPs and HALLMARK gene set (UCell scores)

# Settings ----
## Make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "206_metaprogram_HALLMARK_correlation"

plot_path <- file.path(wd, "plot", analysis_step)
res_path <- file.path(wd, "result", analysis_step)
fp_path <- file.path(wd, "plot", analysis_step, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source(file.path(ws, "common_Rscripts/helpers.R"))

## load data
seu <- readRDS("RDSfiles/seu_203_tum_combined.RDS")

# MP - HALLMARK correlation ----
# parameters ----
sample_col <- "sample"

ucell_cols <- grep("_UCell$", colnames(seu@meta.data), value = TRUE)
mp_cols <- grep("^MP_[0-9]+_UCell$", colnames(seu@meta.data), value = TRUE)
h_cols <- setdiff(ucell_cols, mp_cols)

min_cells_per_sample <- 30
min_mp_sd <- 0.005
min_h_sd <- 0.005

cor_cutoff <- 0.3
fdr_cutoff <- 0.05

# function for one sample ----
cor_mp_h_one_sample <- function(seu_sub) {
  
  sample_id <- unique(seu_sub[[sample_col]][, 1])
  
  if (length(sample_id) != 1) {
    stop("seu_sub should contain exactly one sample")
  }
  
  n_cells <- ncol(seu_sub)
  
  if (n_cells < min_cells_per_sample) {
    return(NULL)
  }
  
  # MP scores from metadata
  mp_df <- seu_sub@meta.data[, mp_cols, drop = FALSE]
  
  # Other UCell scores
  h_df <- seu_sub@meta.data[, h_cols, drop = FALSE]
  
  res <- map_dfr(mp_cols, function(mp) {
    
    mp_score <- mp_df[[mp]]
    
    if (sd(mp_score, na.rm = TRUE) < min_mp_sd) {
      return(NULL)
    }
    
    tmp <- map_dfr(h_cols, function(gene_set) {
      
      h_score <- h_df[, gene_set]
      
      if (sd(h_score, na.rm = TRUE) < min_h_sd) {
        return(NULL)
      }
      
      ct <- suppressWarnings(
        cor.test(
          mp_score,
          h_score,
          method = "spearman",
          exact = FALSE
        )
      )
      
      tibble(
        sample = sample_id,
        n_cells = n_cells,
        MP = mp,
        gene_set = gene_set,
        rho = unname(ct$estimate),
        pval = ct$p.value
      )
    })
    
    # FDR correction within this sample-MP pair
    tmp %>%
      mutate(
        FDR = p.adjust(pval, method = "BH"),
        direction = case_when(
          rho > 0 ~ "positive",
          rho < 0 ~ "negative",
          TRUE ~ "zero"
        )
      )
  })
  
  res
}

# split Seurat object by sample ----
seu_list <- SplitObject(seu, split.by = sample_col)

# run correlation per sample ----
mp_h_cor_sample <- map_dfr(
  seu_list,
  cor_mp_h_one_sample
)

# significant sample-level correlations ----
mp_h_cor_sig <- mp_h_cor_sample %>%
  filter(
    FDR < fdr_cutoff,
    abs(rho) >= cor_cutoff
  )

# summarize recurrence across samples ----
mp_h_cor_summary <- mp_h_cor_sig %>%
  group_by(MP, gene_set, direction) %>%
  summarise(
    n_samples = n_distinct(sample),
    total_cells = sum(n_cells),
    mean_rho = mean(rho, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    min_rho = min(rho, na.rm = TRUE),
    max_rho = max(rho, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(MP, desc(n_samples), desc(abs(mean_rho)))

# recurrent candidates ----
mp_h_cor_recurrent <- mp_h_cor_summary %>%
  filter(n_samples >= 3)

# save ----
write_csv(mp_h_cor_sample, file.path(res_path, "MP_HALLMARK_correlations_by_sample.csv"))
write_csv(mp_h_cor_sig, file.path(res_path, "MP_HALLMARK_correlations_sig_by_sample.csv"))
write_csv(mp_h_cor_summary, file.path(res_path, "MP_HALLMARK_correlations_summary.csv"))
write_csv(mp_h_cor_recurrent, file.path(res_path, "MP_HALLMARK_correlations_recurrent.csv"))
