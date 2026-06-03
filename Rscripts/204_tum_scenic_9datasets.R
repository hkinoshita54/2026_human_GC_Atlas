# Continued from 203_... and pySCENIC output from 9 data sets (060_)
# 2026-05-21, 29
# sample levels correlation between MPs and regulons

# Settings ----
## Make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "204_tum_scenic_9datasets"

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
auc_mtx <- read_tsv(file = "adata/tum_9datasets/auc_mtx.txt") %>% 
  column_to_rownames(var = "...1") %>% 
  as.matrix %>% t()
auc_mtx <- auc_mtx[,Cells(seu)]
seu[["reg"]] <- CreateAssayObject(data = auc_mtx)

## save
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
seu$RNA_snn_res.0.1 <- NULL
seu$RNA_snn_res.0.2 <- NULL
seu$RNA_snn_res.0.3 <- NULL
seu$RNA_snn_res.0.6 <- NULL
seu$annotation2 <- seu$celltype2
seu$celltype2 <- seu$celltype1
seu$celltype1 <- "Epi"
saveRDS(seu, "RDSfiles/seu_203_tum_combined.RDS")

# feature plots ----
DefaultAssay(seu) <- "reg"
add_feat <- c("TEAD1(+)", 
              "ETS2(+)", "FOSL1(+)", "REL(+)",
              "JUN(+)", "FOSB(+)", "ATF3(+)", "FOS(+)", 
              "CEBPG(+)", "IRF7(+)", "IRF7(+)", "PRDM1(+)",
              "FOXM1(+)", "TFDP1(+)", "ILF2(+)", "E2F1(+)", "MYC(+)")
sapply(add_feat, save_fp, seu, fp_path)

# MP - regulon correlation ----
# parameters ----
sample_col <- "sample"

mp_cols <- grep("^MP_[0-9]+_UCell$", colnames(seu@meta.data), value = TRUE)

min_cells_per_sample <- 30
min_mp_sd <- 0.005
min_reg_sd <- 0.005

cor_cutoff <- 0.3
fdr_cutoff <- 0.05

# function for one sample ----
cor_mp_reg_one_sample <- function(seu_sub) {
  
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
  
  # regulon AUC matrix: regulons x cells
  reg_mtx <- LayerData(
    seu_sub,
    assay = "reg",
    layer = "data"
  ) 
  
  # make sure cell order matches
  mp_df <- mp_df[colnames(reg_mtx), , drop = FALSE]
  
  res <- map_dfr(mp_cols, function(mp) {
    
    mp_score <- mp_df[[mp]]
    
    if (sd(mp_score, na.rm = TRUE) < min_mp_sd) {
      return(NULL)
    }
    
    tmp <- map_dfr(rownames(reg_mtx), function(regulon) {
      
      reg_score <- reg_mtx[regulon, ]
      
      if (sd(reg_score, na.rm = TRUE) < min_reg_sd) {
        return(NULL)
      }
      
      ct <- suppressWarnings(
        cor.test(
          mp_score,
          reg_score,
          method = "spearman",
          exact = FALSE
        )
      )
      
      tibble(
        sample = sample_id,
        n_cells = n_cells,
        MP = mp,
        regulon = regulon,
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
mp_reg_cor_sample <- map_dfr(
  seu_list,
  cor_mp_reg_one_sample
)

# significant sample-level correlations ----
mp_reg_cor_sig <- mp_reg_cor_sample %>%
  filter(
    FDR < fdr_cutoff,
    abs(rho) >= cor_cutoff
  )

# summarize recurrence across samples ----
mp_reg_cor_summary <- mp_reg_cor_sig %>%
  group_by(MP, regulon, direction) %>%
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
mp_reg_cor_recurrent <- mp_reg_cor_summary %>%
  filter(n_samples >= 3)

# save ----
write_csv(mp_reg_cor_sample, file.path(res_path, "MP_regulon_correlations_by_sample.csv"))
write_csv(mp_reg_cor_sig, file.path(res_path, "MP_regulon_correlations_sig_by_sample.csv"))
write_csv(mp_reg_cor_summary, file.path(res_path, "MP_regulon_correlations_summary.csv"))
write_csv(mp_reg_cor_recurrent, file.path(res_path, "MP_regulon_correlations_recurrent.csv"))