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
        c("-k", "--k_num"),
        action = "store",
        default = "30",
        type = 'integer',
        help = 'number of nearest neighbors to consider when identifying MNNs.'
    ),
    make_option(
        c("-s", "--sigma"),
        action = "store",
        default = "0.1",
        type = 'integer',
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

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(batchelor))

# args
batch_key <- opt$batch_key
assay_name <- opt$assay_name
k_num <- opt$k_num
n_pcs <- opt$n_pcs
cos_norm <- opt$cos_norm
corrected_assay <- opt$corrected_assay
svd_dim <-  opt$svd_dim

# input file
dataset <- readRDS(opt$input_object)
# batch cell label vector
batch_vector <- as.character(dataset[[opt$batch_key]])
batch_names <- unique(batch_vector)
N_batches <- length(batch_names)
# split object by batches
#batch_list <- lapply(1:N_batches, function(x) {logcounts(dataset[, batch_vector == batch_names[x]])})
batch_list <- lapply(1:N_batches, function(x) {assay(dataset[, batch_vector == batch_names[x]], assay_name)})
# run mnnCorrect
corrected <- do.call('mnnCorrect', c(batch_list, c(k=k_num, sigma=sigma, cos.norm.in=cos_norm, svd.dim=svd_dim)))
# append the corrected expression matrix to object preserving column order
assay(dataset, corrected_assay) <- assay(corrected, "corrected")[, colnames(assay(dataset, assay_name))]
saveRDS(dataset, file = opt$output_object) 
print("mnnCorrect worked!")
