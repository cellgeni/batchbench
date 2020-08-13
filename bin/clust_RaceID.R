#!/usr/bin/env Rscript

# RaceID clustering

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to SCE object input file'
  ),
  make_option(
    c("-f", "--input_features"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Features names subsetted by coefficient of variation'
  ),
  make_option(
    c("-c", "--corrected_assay"),
    action = "store",
    default = "corrected",
    type = 'character',
    help = 'Corrected counts assay name'
  ),
  make_option(
    c("-d", "--dist_metric"),
    action = "store",
    default = "pearson",
    type = 'character',
    help = 'Distance metric'
  ),
  make_option(
    c("-o", "--output_clusters"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to csv output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))


# args
corrected_assay <- opt$corrected_assay
dist_metric <- opt$dist_metric

suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(RaceID))

# read input file
dataset <- readRDS(opt$input_object)
features <- as.character(read.csv(opt$input_features, row.names =1, header = T)$x)
# subset input object by features
dataset <- dataset[features, ]

# 1. Coerce corrected matrix to 'dgCMatrix' class
mat <- as(assay(dataset, corrected_assay), "dgCMatrix")
# 2. Build RaceID object
## 2.1 Remove negative count values (** Should we scale values to zero-minimum instead of removing negative ones?)
mat [mat < 0] <- 0
sc <- SCseq(mat)
# 3. Add features to object (to allow successful execution)
# 3.1 Set expression data (already filtered and normalized) as filtered data (@ndata)
sc@ndata <- sc@expdata
# 3.2 Add total transcript counts
sc@counts <- colSums(as.matrix(sc@expdata))
# 3.3 Add feature names (for distance calculation)
sc@genes <- rownames(sc@ndata)
# 4. Compute distance matrix
sc <- compdist(sc, metric=dist_metric, FSelect = FALSE)
# 5. Cluster. FUNcluster arg must be one of: "kmedoids", "kmeans", "hclust".
# Run all possible clustering types (FUNcluster arg)
sc_kmedoids <- clustexp(sc, clustnr = 30, FUNcluster = "kmedoids", verbose = FALSE)
sc_kmeans <- clustexp(sc, clustnr = 30, FUNcluster = "kmeans", verbose = FALSE)
sc_hclust <- clustexp(sc, clustnr = 30, FUNcluster = "hclust", verbose = FALSE)
# 6. Save cluster annotation
write.csv(data.frame("kmedoids.RaceID" = sc_kmedoids@cluster[["kpart"]]), file = paste0("kmedoids_", opt$output_clusters))
write.csv(data.frame("kmeans.RaceID" = sc_kmeans@cluster[["kpart"]]), file = paste0("kmeans_", opt$output_clusters))
write.csv(data.frame("hclust.RaceID" = sc_hclust@cluster[["kpart"]]), file = paste0("hclust_", opt$output_clusters))

print("RaceID clustering annotations saved")
