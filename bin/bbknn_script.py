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

def main(args):
    '''This function reads an h5ad file, runs PCA, runs BBKNN and saves the corrected object'''

    #read input h5ad
    dataset = sc.read(args.input)
    print("File read!")
    #compute PCA
    sc.tl.pca(dataset, n_comps = 10)
    #run bbknn
    dataset_bbknn = bbknn.bbknn(dataset,batch_key='Batch', neighbors_within_batch= 10, n_pcs=10, save_knn=True, copy=True)
    print("BBKNN done!")

    #save corrected object to the specified path in the output
    dataset_bbknn.write(args.output)


#this is the main entry point to the compiler to go through when reading the script
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument('input',
                        help='input file to be batch corrected')

    parser.add_argument('output',
                        help='output file batch corrected')

    args = parser.parse_args()

    main(args)
