# 2026-05-18
# NMF based metaprogram discovery as in Gavish et al. nature 2023, continued from 200.1_...R
# https://github.com/tiroshlab/3ca/tree/main/ITH_hallmarks

# Settings ----
## make directories
wd <- getwd()
ws <- "/Users/hiroto/WORKSPACE"

analysis_step <- "201.1_metaprogram"

plot_path <- file.path(wd, "plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path(wd, "result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)

# load data ----
Genes_nmf_w_basis <- readRDS("RDSfiles/Genes_nmf_w_basis_log2CPM10_centered_rank4_9_nrun10_min10cells.RDS")
Genes_nmf_w_basis_2 <- readRDS("RDSfiles/Genes_nmf_w_basis_log2CPM10_centered_rank4_9_nrun10_min10cells_additional.RDS")
Genes_nmf_w_basis <- c(Genes_nmf_w_basis, Genes_nmf_w_basis_2)

# generate metaprogram ----

## 0. Helper functions ----
get_patient_from_program <- function(x) {sub("_[0-9]+_[0-9]+$", "", x)}

overlap_mat <- function(mat) {
  apply(mat, 2, function(x) {
    apply(mat, 2, function(y) length(intersect(x, y)))
  })
}

robust_nmf_programs <- function(nmf_programs,
                                intra_min = 35,
                                intra_max = 10,
                                inter_filter = TRUE,
                                inter_min = 10) {
  ## nmf_programs: list of matrices, one per patient/sample
  ## each matrix: top_genes x programs, entries are gene symbols
  
  ## 1) Robust within patient: program must overlap strongly with another program from the same patient.
  intra_intersect <- lapply(nmf_programs, overlap_mat)
  
  intra_intersect_max <- lapply(intra_intersect, function(x) {
    apply(x, 2, function(y) sort(y, decreasing = TRUE)[2])
  })
  
  nmf_sel <- lapply(names(nmf_programs), function(pid) {
    z <- nmf_programs[[pid]]
    keep <- intra_intersect_max[[pid]] >= intra_min
    z[, keep, drop = FALSE]
  })
  names(nmf_sel) <- names(nmf_programs)
  
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  if (is.null(nmf_sel_unlist) || ncol(nmf_sel_unlist) == 0) {
    stop("No NMF programs passed intra-sample robustness filtering.")
  }
  
  ## 2) Inter-patient similarity and non-redundancy within patient
  inter_intersect <- overlap_mat(nmf_sel_unlist)
  
  final_filter <- character()
  
  for (pid in names(nmf_sel)) {
    pid_cols <- colnames(nmf_sel[[pid]])
    other_cols <- setdiff(colnames(nmf_sel_unlist), pid_cols)
    
    if (length(pid_cols) == 0 || length(other_cols) == 0) next
    
    a <- inter_intersect[other_cols, pid_cols, drop = FALSE]
    b <- sort(apply(a, 2, max), decreasing = TRUE)
    
    if (inter_filter) {
      b <- b[b >= inter_min]
    }
    
    if (length(b) == 0) next
    
    ## Iteratively keep programs with high inter-patient similarity, while avoiding redundant programs from the same patient.
    chosen <- names(b)[1]
    
    if (length(b) > 1) {
      for (j in 2:length(b)) {
        candidate <- names(b)[j]
        if (max(inter_intersect[chosen, candidate]) <= intra_max) {
          chosen <- c(chosen, candidate)
        }
      }
    }
    
    final_filter <- c(final_filter, chosen)
  }
  
  unique(final_filter)
}

get_nmf_scores_for_genes <- function(genes, programs, Genes_nmf_w_basis) {
  ## Used for tie-breaking around the 50th gene.
  scores <- c()
  
  for (prog in programs) {
    pid <- get_patient_from_program(prog)
    mat <- Genes_nmf_w_basis[[pid]]
    
    if (is.null(mat) || !prog %in% colnames(mat)) next
    
    gene_match <- match(genes, toupper(rownames(mat)))
    ok <- !is.na(gene_match)
    
    q <- mat[gene_match[ok], prog]
    names(q) <- genes[ok]
    scores <- c(scores, q)
  }
  
  scores
}

update_MP_genes <- function(NMF_history,
                            Curr_cluster,
                            Genes_nmf_w_basis,
                            top_n = 50) {
  tab <- sort(table(NMF_history), decreasing = TRUE)
  
  if (length(tab) <= top_n) {
    return(names(tab))
  }
  
  border_freq <- tab[top_n]
  genes_above <- names(tab[tab > border_freq])
  genes_at_border <- names(tab[tab == border_freq])
  
  n_needed <- top_n - length(genes_above)
  
  if (length(genes_at_border) > n_needed) {
    scores <- get_nmf_scores_for_genes(
      genes = genes_at_border,
      programs = Curr_cluster,
      Genes_nmf_w_basis = Genes_nmf_w_basis
    )
    
    scores <- sort(scores, decreasing = TRUE)
    scores <- scores[!duplicated(names(scores))]
    
    genes_border_sorted <- names(scores)
    
    ## In case some tied genes were not recovered in scores
    genes_border_sorted <- c(
      genes_border_sorted,
      setdiff(genes_at_border, genes_border_sorted)
    )
    
    c(genes_above, genes_border_sorted)[1:top_n]
  } else {
    names(tab)[1:top_n]
  }
}

## 1. Safety checks ----
stopifnot(is.list(Genes_nmf_w_basis))
stopifnot(!is.null(names(Genes_nmf_w_basis)))

## Make sure list names match patient IDs in column names
all_colnames <- unlist(lapply(Genes_nmf_w_basis, colnames))
all_pids_from_cols <- unique(get_patient_from_program(all_colnames))

missing_pids <- setdiff(all_pids_from_cols, names(Genes_nmf_w_basis))
if (length(missing_pids) > 0) {
  stop("These patient IDs are present in column names but not list names: ",
       paste(missing_pids, collapse = ", "))
}

## 2. Convert each NMF program to top 50 genes
top_n_genes <- 50

nmf_programs <- lapply(Genes_nmf_w_basis, function(x) {
  apply(x, 2, function(y) {
    toupper(names(sort(y, decreasing = TRUE))[seq_len(top_n_genes)])
  })
})

## Ensure matrices stay matrices
nmf_programs <- lapply(nmf_programs, as.matrix)

## 3. Select robust NMF programs
intra_min_parameter <- 35
intra_max_parameter <- 10
inter_min_parameter <- 10

nmf_filter <- robust_nmf_programs(
  nmf_programs,
  intra_min = intra_min_parameter,
  intra_max = intra_max_parameter,
  inter_filter = TRUE,
  inter_min = inter_min_parameter
)

message("Selected robust NMF programs: ", length(nmf_filter))

nmf_programs_filt_list <- lapply(nmf_programs, function(x) {
  x[, colnames(x) %in% nmf_filter, drop = FALSE]
})

nmf_programs_filt <- do.call(cbind, nmf_programs_filt_list)

if (ncol(nmf_programs_filt) == 0) {
  stop("No programs remained after robust filtering.")
}

## 4. Pairwise top-50 overlap
nmf_intersect <- overlap_mat(nmf_programs_filt)

nmf_intersect_hc <- hclust(as.dist(top_n_genes - nmf_intersect), method = "average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))

ord <- order.dendrogram(nmf_intersect_hc)
nmf_intersect <- nmf_intersect[ord, ord]

nmf_intersect_original <- nmf_intersect
nmf_programs_current <- nmf_programs_filt[, colnames(nmf_intersect), drop = FALSE]

## 5. Iteratively generate MPs
Min_intersect_initial <- 10
Min_intersect_cluster <- 10
Min_group_size <- 5

Cluster_list <- list()
MP_list <- list()

k <- 1

Sorted_intersection <- sort(
  apply(nmf_intersect, 2, function(x) length(which(x >= Min_intersect_initial)) - 1),
  decreasing = TRUE
)

while (length(Sorted_intersection) > 0 &&
       Sorted_intersection[1] > Min_group_size &&
       ncol(nmf_programs_current) > 1) {
  
  seed_prog <- names(Sorted_intersection)[1]
  Curr_cluster <- seed_prog
  
  Genes_MP <- nmf_programs_current[, seed_prog]
  NMF_history <- Genes_MP
  
  ## Remove seed program
  nmf_programs_current <- nmf_programs_current[
    ,
    setdiff(colnames(nmf_programs_current), seed_prog),
    drop = FALSE
  ]
  
  if (ncol(nmf_programs_current) == 0) break
  
  Intersection_with_Genes_MP <- sort(
    apply(nmf_programs_current, 2, function(x) length(intersect(Genes_MP, x))),
    decreasing = TRUE
  )
  
  while (length(Intersection_with_Genes_MP) > 0 &&
         Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {
    
    next_prog <- names(Intersection_with_Genes_MP)[1]
    
    Curr_cluster <- c(Curr_cluster, next_prog)
    
    NMF_history <- c(NMF_history, nmf_programs_current[, next_prog])
    
    Genes_MP <- update_MP_genes(
      NMF_history = NMF_history,
      Curr_cluster = Curr_cluster,
      Genes_nmf_w_basis = Genes_nmf_w_basis,
      top_n = top_n_genes
    )
    
    ## Remove added program
    nmf_programs_current <- nmf_programs_current[
      ,
      setdiff(colnames(nmf_programs_current), next_prog),
      drop = FALSE
    ]
    
    if (ncol(nmf_programs_current) == 0) break
    
    Intersection_with_Genes_MP <- sort(
      apply(nmf_programs_current, 2, function(x) length(intersect(Genes_MP, x))),
      decreasing = TRUE
    )
  }
  
  Cluster_list[[paste0("Cluster_", k)]] <- Curr_cluster
  MP_list[[paste0("MP_", k)]] <- Genes_MP
  
  message(
    paste0(
      "MP_", k, ": ",
      length(Curr_cluster), " programs from ",
      length(unique(get_patient_from_program(Curr_cluster))), " patients"
    )
  )
  
  ## Remove current cluster from intersection matrix
  remaining <- setdiff(colnames(nmf_intersect), Curr_cluster)
  
  if (length(remaining) <= 1) break
  
  nmf_intersect <- nmf_intersect[remaining, remaining, drop = FALSE]
  
  Sorted_intersection <- sort(
    apply(nmf_intersect, 2, function(x) length(which(x >= Min_intersect_initial)) - 1),
    decreasing = TRUE
  )
  
  k <- k + 1
}

## 6. Outputs
MP_list <- lapply(MP_list, unique)

MP_mat <- do.call(cbind, MP_list)

Cluster_summary <- tibble(
  MP = names(Cluster_list),
  n_programs = lengths(Cluster_list),
  n_patients = sapply(Cluster_list, function(x) length(unique(get_patient_from_program(x)))),
  programs = sapply(Cluster_list, paste, collapse = ";")
)

MP_genes_long <- enframe(MP_list, name = "MP", value = "gene") %>%
  tidyr::unnest(gene) %>%
  group_by(MP) %>%
  mutate(rank = row_number()) %>%
  ungroup()

## Useful objects:
## MP_list         : list; each MP is a vector of top 50 genes
## MP_mat          : 50 x number_of_MPs matrix
## Cluster_list    : NMF programs contributing to each MP
## Cluster_summary : summary table
## MP_genes_long   : long-format MP gene table
## nmf_intersect_original : pairwise top-50 overlap matrix before clustering

saveRDS(
  list(
    MP_list = MP_list,
    MP_mat = MP_mat,
    Cluster_list = Cluster_list,
    Cluster_summary = Cluster_summary,
    MP_genes_long = MP_genes_long,
    nmf_filter = nmf_filter,
    nmf_intersect_original = nmf_intersect_original
  ),
  file = "RDSfiles/meta_program_results_2.rds"
)


