# 2026-05-19
# Inspect the MPs from 201.1_...R

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "202_inspect_metaprograms"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)

## load data
mp_res <- readRDS("RDSfiles/meta_program_results_2.rds")
Cluster_list <- mp_res$Cluster_list
seu <- readRDS("RDSfiles/seu_010_tum.RDS")
seu2 <- readRDS("RDSfiles/seu_010.1_tum2.RDS")
seu <- merge(x = seu, y = seu2)
seu <- JoinLayers(seu)

## helper
get_sample_from_program <- function(x) {
  str_remove(x, "_[0-9]+_[0-9]+$")
}

# Sample metadata ----
sample_meta <- seu@meta.data %>% distinct(sample, dataset)

# Convert Cluster_list to long table ----
mp_program_long <- enframe(Cluster_list, name = "MP", value = "program") %>%
  unnest(program) %>%
  mutate(
    program = as.character(program),
    sample = get_sample_from_program(program)
  ) %>%
  left_join(sample_meta, by = "sample")

# Check unmatched samples ----
unmatched <- mp_program_long %>%
  filter(is.na(dataset)) %>%
  distinct(sample)

if (nrow(unmatched) > 0) {
  warning("Some program samples were not found in Seurat metadata:")
  print(unmatched)
}

# Summary: MP x dataset ----
mp_dataset_summary <- mp_program_long %>%
  group_by(MP, dataset) %>%
  summarise(
    n_programs = n(),
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  group_by(MP) %>%
  mutate(
    total_programs = sum(n_programs),
    frac_programs = n_programs / total_programs
  ) %>%
  ungroup() %>%
  arrange(MP, desc(n_programs))

# Wide table: absolute program counts ----
mp_dataset_n_programs_wide <- mp_dataset_summary %>%
  select(MP, dataset, n_programs) %>%
  pivot_wider(
    names_from = dataset,
    values_from = n_programs,
    values_fill = 0
  )

# Wide table: fraction of programs ----
mp_dataset_frac_wide <- mp_dataset_summary %>%
  select(MP, dataset, frac_programs) %>%
  pivot_wider(
    names_from = dataset,
    values_from = frac_programs,
    values_fill = 0
  )

# Dataset diversity score per MP ----
mp_dataset_diversity <- mp_dataset_summary %>%
  group_by(MP) %>%
  summarise(
    n_datasets = n_distinct(dataset),
    max_dataset_fraction = max(frac_programs),
    dominant_dataset = dataset[which.max(frac_programs)],
    shannon_entropy = -sum(frac_programs * log(frac_programs)),
    total_programs = first(total_programs),
    .groups = "drop"
  ) %>%
  mutate(
    dataset_specific_flag = case_when(
      max_dataset_fraction >= 0.90 ~ "dataset_specific",
      max_dataset_fraction >= 0.70 ~ "dataset_enriched",
      TRUE ~ "multi_dataset"
    )
  ) %>%
  arrange(desc(max_dataset_fraction))

# Save outputs ----
write_csv(mp_program_long, file.path(res_path, "MP_program_long_with_dataset.csv"))
write_csv(mp_dataset_summary, file.path(res_path, "MP_dataset_summary_long.csv"))
write_csv(mp_dataset_n_programs_wide, file.path(res_path, "MP_dataset_n_programs_wide.csv"))
write_csv(mp_dataset_frac_wide, file.path(res_path, "MP_dataset_fraction_wide.csv"))
write_csv(mp_dataset_diversity, file.path(res_path, "MP_dataset_diversity.csv"))

# Quick inspection ----
mp_dataset_diversity
mp_dataset_summary %>% filter(MP %in% c("Cluster_12", "Cluster_15"))

# Plot: fraction of programs by dataset ----
ggplot(mp_dataset_summary, aes(x = MP, y = frac_programs, fill = dataset)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Fraction of contributing NMF programs") +
  xlab(NULL)

ggsave(
  file.path(plot_path, "MP_dataset_fraction_stacked_barplot.pdf"),
  width = 8,
  height = 4
)

# Plot: absolute number of programs by dataset ----
ggplot(mp_dataset_summary, aes(x = MP, y = n_programs, fill = dataset)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of contributing NMF programs") +
  xlab(NULL)

ggsave(
  file.path(plot_path, "MP_dataset_n_programs_stacked_barplot.pdf"),
  width = 8,
  height = 4
)
