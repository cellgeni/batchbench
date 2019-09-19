#!/usr/bin/env python3

#modules
import numpy as np
import pandas as pd
import scanpy.api as sc 

import anndata
import os
import scipy
import sys
import argparse

def save_h5ad(dataset):
    dataset.write(args.output)

def remove_elements(dataset, filter_vec, threshold, info):
    """Function to remove those batches or cell types with less than threshold cells"""
    
    removal_names = filter_vec.value_counts()[filter_vec.value_counts() <= threshold].index.tolist()
    
    if len(removal_names) == 0:
        print("There are no", info, "with less than", threshold, "cells:")   
    if len(removal_names) != 0:
        print("There are", len(removal_names), info, "with less than", threshold, "cells:", removal_names, "removed!")  
   
        for i in range(0, len(removal_names)):
            dataset = dataset[dataset.obs["cell_type1"] != removal_names[i]]
            dataset = dataset[dataset.obs["Batch"] != removal_names[i]]
            dataset = anndata.AnnData.copy(self = dataset)

    return dataset

def do_the_filtering(dataset, **kwargs):
    #remove thos cells with no cell type annotation
    dataset = dataset[dataset.obs['cell_type1'] != "NA"]
    #remove those batches/cell types with less than threshold cells <------FUNCITON!
    print("Dimensions pre-filtering:", dataset.shape)
    dataset = remove_elements(dataset, filter_vec =  kwargs["batch_vector"], threshold = args.batch_threshold , info = "Batch/es")
    dataset = remove_elements(dataset, filter_vec = kwargs["cell_type_vector"], threshold = args.cell_type_threshold, info = "Cell type/s")
    print("Dimensions post-filtering:", dataset.shape)
    
    save_h5ad(dataset)

def read_h5ad(x):    
    #read input h5ad
    dataset = sc.read(x)
    
    kwargs = {}
    kwargs["batch_vector"] = dataset.obs["Batch"]
    kwargs["cell_type_vector"] = dataset.obs["cell_type1"]

    do_the_filtering(dataset, **kwargs)

if __name__== "__main__":
    
    parser = argparse.ArgumentParser(description='Input/Output files and other parameters')

    parser.add_argument("--input", dest='input',
                        help ='h5ad object where some QC is going to be applied')
    
    parser.add_argument("--batch_threshold", dest='batch_threshold', type = int,
                        help ='minimum number of cells required for a BATCH to be present')
        
    parser.add_argument("--cell_type_threshold", dest='cell_type_threshold', type = int,
                        help ='minimum number of cells required for a CELL TYPE to be present')
    
    parser.add_argument('--output', dest='output',
                        help='QCded h5ad object')

    args = parser.parse_args()

    read_h5ad(args.input)
