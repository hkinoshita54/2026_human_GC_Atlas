# to check common features across data sets
# 2026-03-28

library(Seurat)
seu <- readRDS("RDSfiles/seu_010_tum.RDS")

# check dataset names
table(seu$dataset)

# counts matrix from RNA assay
cts <- GetAssayData(seu, assay = "RNA", layer = "counts")

# split cell names by dataset
cells.by.dataset <- split(colnames(seu), seu$dataset)

# genes detected in each dataset (nonzero in at least one cell)
features.by.dataset <- lapply(cells.by.dataset, function(cells) {
  rownames(cts)[Matrix::rowSums(cts[, cells, drop = FALSE] > 0) > 0]
})

# number of detected genes in each dataset
sapply(features.by.dataset, length)

# pairwaise overlaps 
datasets <- names(features.by.dataset)

overlap.mat <- outer(
  datasets, datasets,
  Vectorize(function(x, y) length(intersect(features.by.dataset[[x]], features.by.dataset[[y]])))
)

dimnames(overlap.mat) <- list(datasets, datasets)
overlap.mat

# intersection of all
common.features <- Reduce(intersect, features.by.dataset)
length(common.features)
head(common.features)

# unique features
unique.features <- lapply(datasets, function(d) {
  others <- setdiff(datasets, d)
  setdiff(features.by.dataset[[d]], Reduce(union, features.by.dataset[others]))
})
names(unique.features) <- datasets

sapply(unique.features, length)