#!/usr/bin/env python3

import argparse
import bbknn
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scipy

# save bbknn output
def save_corrected_object(dataset):
    dataset.write(args.output)
  
def main(dataset):
    '''This function runs PCA, runs BBKNN and saves the corrected object'''
    # compute PCA
    sc.tl.pca(dataset, n_comps = 15)
    # run bbknn
    dataset_bbknn = bbknn.bbknn(dataset, 
				batch_key=args.batch_key, 
				neighbors_within_batch=args.n_neighbours, 
				n_pcs=args.n_pcs, 
				copy=True)
    print("BBKNN done!")

    save_corrected_object(dataset_bbknn)

# read input object
def read_h5ad(dataset):
    return sc.read(dataset)
    print("File read!")
 
# args
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='Input h5ad object.')
    parser.add_argument("--batch_key", dest='batch_key',
                        help ='Cell key defining Batch')
    parser.add_argument("--n_pcs", dest='n_pcs',
                        help ='Number of PCs for PCA prior to graph correction.')
    parser.add_argument("--n_neighbours", dest='n_neighbours',
                        help ='Number of nearest neighbours per batch to perform graph correction.')
    parser.add_argument('--output', dest='output',
                        help='AnnData object with corrected graph in object.uns["neighbours"].') 

    args = parser.parse_args()

    main(read_h5ad(args.input))
