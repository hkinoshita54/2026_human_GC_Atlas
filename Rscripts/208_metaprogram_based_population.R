# Continued from 204_... 
# 2026-05-29, 31
# Define cell populations based on metaprograms

# Settings ----
## Make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "208_metaprogram_based_populations"

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

# Define MP based populations ----
# ## MP4
# seu$MP4_q <- ntile(seu$MP_4_UCell, 4) # divide scores into quantiles
# seu$MP4_state <- case_when(
#   seu$MP4_q == 4 ~ "MP4_4th",
#   seu$MP4_q %in% c(1,2) ~ "MP4_1st2nd",
#   TRUE ~ "MP4_3rd"
# )
# 
# ## MP5
# seu$MP5_q <- ntile(seu$MP_5_UCell, 4) # divide scores into quantiles
# seu$MP5_state <- case_when(
#   seu$MP5_q == 4 ~ "MP5_4th",
#   seu$MP5_q %in% c(1,2) ~ "MP5_1st2nd",
#   TRUE ~ "MP5_3rd"
# )
# 
# ## MP14
# seu$MP14_q <- ntile(seu$MP_14_UCell, 4) # divide scores into quantiles
# seu$MP14_state <- case_when(
#   seu$MP14_q == 4 ~ "MP14_4th",
#   seu$MP14_q %in% c(1,2) ~ "MP14_1st2nd",
#   TRUE ~ "MP14_3rd"
# )

# scoring with sequencing depth correction ----
seu$log_nFeature <- log10(seu$nFeature_RNA + 1)

## MP4 with sequencing depth correction
fit <- lm(MP_4_UCell ~ log_nFeature + dataset, data = seu@meta.data)
seu$MP4_resid <- residuals(fit)
seu$MP4_resid_q <- ntile(seu$MP4_resid, 4)
seu$MP4_resid_q <- case_when(
  seu$MP4_resid_q == 4 ~ "MP4_4",
  seu$MP4_resid_q %in% c(1,2) ~ "MP4_1_2",
  TRUE ~ "MP4_3"
)

## MP5 with sequencing depth correction
fit <- lm(MP_5_UCell ~ log_nFeature + dataset, data = seu@meta.data)
seu$MP5_resid <- residuals(fit)
seu$MP5_resid_q <- ntile(seu$MP5_resid, 4)
seu$MP5_resid_q <- case_when(
  seu$MP5_resid_q == 4 ~ "MP5_4",
  seu$MP5_resid_q %in% c(1,2) ~ "MP5_1_2",
  TRUE ~ "MP5_3"
)

## MP14 with sequencing depth correction
fit <- lm(MP_14_UCell ~ log_nFeature + dataset, data = seu@meta.data)
seu$MP14_resid <- residuals(fit)
seu$MP14_resid_q <- ntile(seu$MP14_resid, 4)
seu$MP14_resid_q <- case_when(
  seu$MP14_resid_q == 4 ~ "MP14_4",
  seu$MP14_resid_q %in% c(1,2) ~ "MP14_1_2",
  TRUE ~ "MP14_3"
)

DimPlot(seu, group.by = "MP4_resid_q", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "MP4resid_q")
ggsave("MP4resid_q.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "MP5_resid_q", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "MP5resid_q")
ggsave("MP5resid_q.pdf", path = plot_path, width = 35, height = 30, units = "mm")

DimPlot(seu, group.by = "MP14_resid_q", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "MP14resid_q")
ggsave("MP14resid_q.pdf", path = plot_path, width = 35, height = 30, units = "mm")

saveRDS(seu, "RDSfiles/seu_203_tum_combined.RDS")

# Sample level analyses ----
dat_all <- seu@meta.data %>%
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
    n_cells = n(),
    
    n_MP4_high = sum(MP4_resid_q == "MP4_4"),
    prop_MP4_high = n_MP4_high / n_cells,
    mean_MP4_score = mean(MP4_resid),
    
    n_MP5_high = sum(MP5_resid_q == "MP5_4"),
    prop_MP5_high = n_MP5_high / n_cells,
    mean_MP5_score = mean(MP5_resid),
    
    n_MP14_high = sum(MP14_resid_q == "MP14_4"),
    prop_MP14_high = n_MP14_high / n_cells,
    mean_MP14_score = mean(MP14_resid),
  ) 

## Diffuse vs Intestinal
dat <- dat_all %>% filter(Lauren %in% c("Diffuse", "Intestinal"))
dat$Lauren <- factor(dat$Lauren, levels = c("Intestinal", "Diffuse"))

fit <- lm(as.formula(paste0("prop_MP4_high", " ~ Lauren + dataset")), data = dat)
summary(fit)

fit <- lm(as.formula(paste0("prop_MP5_high", " ~ Lauren + dataset")), data = dat)
summary(fit)

fit <- lm(as.formula(paste0("prop_MP14_high", " ~ Lauren + dataset")), data = dat)
summary(fit)

# TCGA
dat <- dat_all %>% filter(!is.na(TCGA))
dat$TCGA <- factor(dat$TCGA, levels = c("CIN", "GS", "MSI", "EBV"))

fit <- lm(as.formula(paste0("prop_MP4_high", " ~ TCGA + dataset")), data = dat)
summary(fit)

fit <- lm(as.formula(paste0("prop_MP5_high", " ~ TCGA + dataset")), data = dat)
summary(fit)

fit <- lm(as.formula(paste0("prop_MP14_high", " ~ TCGA + dataset")), data = dat)
summary(fit)

# ACRG
dat <- dat_all %>% filter(!is.na(ACRG))
dat$ACRG <- factor(dat$ACRG, levels = c("MSS/TP53-", "EMT", "MSI", "MSS/TP53+"))

fit <- lm(as.formula(paste0("prop_MP4_high", " ~ ACRG + dataset")), data = dat)
summary(fit)

fit <- lm(as.formula(paste0("prop_MP5_high", " ~ ACRG + dataset")), data = dat)
summary(fit)

fit <- lm(as.formula(paste0("prop_MP14_high", " ~ ACRG + dataset")), data = dat)
summary(fit)
