#!/usr/bin/env python3
  
#modules
import bbknn
import numpy as np
import pandas as pd
import scanpy.api as sc
import anndata
import os
import scipy
import sys

import argparse
def save_corrected_object(dataset):
    #save corrected object to the specified path in the output
    dataset.write(args.output)
  
def main(dataset):
    '''This function runs PCA, runs BBKNN and saves the corrected object'''
    #compute PCA
    sc.tl.pca(dataset, n_comps = 10)
    #run bbknn
    dataset_bbknn = bbknn.bbknn(dataset,batch_key='Batch', neighbors_within_batch= 10, n_pcs=10, save_knn=True, copy=True)
    print("BBKNN done!")

    save_corrected_object(dataset_bbknn)

def read_h5ad(dataset):
    #read input h5ad
    return sc.read(dataset)
    print("File read!")
 
    
#this is the main entry point to the compiler to go through when reading the script
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='h5ad object over which BBKNN is going to be ran')

    parser.add_argument('--output', dest='output',
                        help='BBKNN corrected object stored as h5ad.BBKNN corrects the kNN graph, so corrected embedding is stored in object.uns[neighbours]')

    args = parser.parse_args()

    main(read_h5ad(args.input))
