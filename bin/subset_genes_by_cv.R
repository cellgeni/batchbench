#!/usr/bin/env Rscript

# Rank genes by its coefficient of variation and subset. 
# This is to consider different amount of genes in the clustering analysis.

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
    c("-a", "--assay_name"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Counts assay to add to the h5ad object'
  ),
  make_option(
    c("-p", "--prop_genes"),
    action = "store",
    default = 1.0,
    type = 'numeric',
    help = 'Proportion of genes to subset by their coefficient of variation.'
  ),
  make_option(
    c("-o", "--output_genes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to csv output file with gene names ranked by '
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# FUNCTIONS #
coef_var <- function(gene_exp){
  cv <- sd(gene_exp)/mean(gene_exp)*100
  cv
}

# args
assay_name <- opt$assay_name
prop_genes <- opt$prop_genes

suppressPackageStartupMessages(require(SingleCellExperiment))
# read input file
dataset <- readRDS(opt$input_object)
# compute coefficient of variation
gene_coef_var <- apply(assay(dataset, assay_name), 1, coef_var)
# rank by CV
rank_gene_coef_var <- sort(gene_coef_var, decreasing = T)
# absolute number of genes to retain
n_genes <- round(nrow(dataset)* prop_genes, 0)
# subset genes by coef variation
sub_genes_cv <- names(rank_gene_coef_var[1: n_genes])
# save names of genes subsetted by their coefficient of variation
write.csv(sub_genes_cv, opt$output_genes)
