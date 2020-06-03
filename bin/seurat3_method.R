#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list = list(
    make_option(
        c("-i", "--input_object"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to rds input file' 
    ),
    make_option(
        c("-b", "--batch_key"),
        action = "store",
        default = "Batch",
        type = 'character',
        help = 'Minimum number of cells for a gene to be expressed in.'
    ),
    make_option(
        c("-h", "--hvg_method"),
        action = "store",
        default = "dispersion",
        type = 'character',
        help = 'Select method to FindVariableFeatures: one of "vst", "mean.var.plot" or "dispersion"' 
    ),
    make_option(
        c("-f", "--n_features"),
        action = "store",
        default = "2000",
        type = 'integer',
        help = 'N of variable features to find in FindVariableFeatures' 
    ),
    make_option(
        c("-r", "--anchors"),
        action = "store",
        default = "30",
        type = 'integer',
        help = ' Number of anchors to find to perform integration'
    ),
    make_option(
        c("-s", "--sigma"),
        action = "store",
        default = "0.1",
        type = 'numeric',
        help = 'numeric scalar specifying the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell'
    ),
    make_option(
        c("-v", "--svd_dim"),
        action = "store",
        default = 2,
        type = 'integer',
        help = 'number of dimensions to use for summarizing biological substructure within each batch.'
    ),
    make_option(
        c("-a", "--assay_name"),
        action = "store",
        default = "logcounts",
        type = 'character',
        help = 'Counts assay to add to the h5ad object'
    ),
    make_option(
        c("-c", "--corrected_assay"),
        action = "store",
        default = "corrected",
        type = 'character',
        help = 'Corrected counts assay name'
    ),
    make_option(
        c("-n", "--cos_norm"),
        action = "store",
        default = TRUE,
        type = 'logical',
        help = ' should cosine normalization be performed on the input data prior to PCA.' 
    ),
    make_option(
        c("-o", "--output_object"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to rds output file'
    )
)

print("I")
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Seurat))
# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay 
batch_key <- opt$batch_key
hvg_method <- opt$hvg_method
n_features <- opt$n_features
n_anchors <- opt$n_anchors

# input dataset
dataset <- readRDS(opt$input_object)
# convert to seurat object
dataset <- as.Seurat(dataset, counts = assay_name, data = assay_name)
batch_vector <- as.character(dataset[[batch_key]])
batch_names  <- names(table(batch_vector))
N_batches <- length(batch_names)

# split seurat object into batches
batch_list <- lapply(1:N_batches, function(x) {abc <- dataset[, dataset[[batch_key]] == batch_names[x]]})
# Normalize and find HVG
for (i in 1:N_batches) {
  batch_list[[i]] <- NormalizeData(object = batch_list[[i]], 
						verbose = FALSE)
  batch_list[[i]] <- FindVariableFeatures(object = batch_list[[i]], 
                                          	selection.method = hvg_method, 
						nfeatures = n_features, 
						verbose = F)
	}

# prevent small datasets from not having enough neighbors (k) to use when filtering anchors 
if(any(sapply(batch_list, ncol)) < 200) {
  k_filter <- (min(sapply(batch_list, ncol)))
  }else{k_filter =200}
# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:n_anchors, k.filter = k_filter)
# Integrate subsets
integrated <- IntegrateData( new.assay.name = corrected_assay, anchorset = anchors, dims = 1:n_anchors)
# save seurat3 corrected object
saveRDS(integrated, opt$output_object)
