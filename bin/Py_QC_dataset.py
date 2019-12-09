#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc 
import anndata
import argparse

def save_h5ad(dataset):
    dataset.write(args.output)

def print_filtering(dataset, filter_vec,  threshold, meta_name):
    """Function to select the filtering_names(names of those batches or cell types with less proportion of cells than threshold),
    and print an informative table with: batches/cell types, absolute_n_cells, relative_n_cells, Exluded or not.
    """
    cell_count = filter_vec.value_counts(ascending=False)
    
    print("**", meta_name ,  "containing less than:",  str(threshold),  "of total cells are removed" +"\n" + "**", meta_name, "filtered based on our threshold")
    
    #dataframe informing about the filtering about to be done
    exclude_df = pd.DataFrame({meta_name: cell_count.index.to_list(), 'n_cells': cell_count.values, 
              '%_cells': cell_count.values/dataset.n_obs, 'Excluded_?': cell_count.values/dataset.n_obs < threshold})
    print(exclude_df)
    
    removal_names = exclude_df[meta_name][exclude_df["Excluded_?"] == True].tolist()
    return removal_names 

def pre_processing(dataset, min_genes, min_cells):
    """QC pre processing funciton to remove: genes expressed in less than min_cells, 
    and cells with less than min_genes expressed"""
    sc.pp.filter_cells(dataset, min_genes = min_genes)
    print("** CELLS with less than", str(min_genes), "genes expressed filtered out.")
    sc.pp.filter_genes(dataset, min_cells = min_cells)
    print("** GENES expressed in less than", str(min_cells), "cells filtered out.")
    
    print("** Dimensions after preprocessing:", dataset.shape)
    return dataset

def remove_elements(dataset, filter_vec, threshold, info_meta, meta_name):
    """Function to remove those batches or cell types with less than threshold % cells"""
    #remove thos cells with no cell type annotation
    dataset = dataset[dataset.obs[info_meta] != "NA"]
    print("** Cells of", meta_name, " annotated as NA removed")
    
    removal_names = print_filtering(dataset, filter_vec, threshold, meta_name)

    for i in range(0, len(removal_names)):
        dataset = dataset[dataset.obs["cell_type1"] != removal_names[i]]
        dataset = dataset[dataset.obs["Batch"] != removal_names[i]]
        dataset = anndata.AnnData.copy(self = dataset)
    return dataset
    
def do_the_filtering(dataset, **kwargs):
    #pre processing(remove cells with les than X genes, and genes expressed in less than Y cells)
    print("** Dimensions pre-filtering:", dataset.shape)
    dataset = pre_processing(dataset, args.min_genes, args.min_cells)
    
    #remove those batches/cell types with less than threshold cells <------FUNCITON!
    dataset = remove_elements(dataset, filter_vec =  dataset.obs["Batch"], threshold = args.bt_thres , info_meta = "Batch", meta_name = "Batch/es")
    dataset = remove_elements(dataset, filter_vec = dataset.obs["cell_type1"], threshold = args.ct_thres, info_meta = "cell_type1", meta_name = "Cell type/s")
    print("** Dimensions post-filtering:", dataset.shape)
    
    save_h5ad(dataset)

def read_h5ad(x):    
    #read input h5ad
    dataset = sc.read(x)
    kwargs = {}
    #kwargs["batch_vector"] = dataset.obs["Batch"]
    #kwargs["cell_type_vector"] = dataset.obs["cell_type1"]

    do_the_filtering(dataset, **kwargs)

if __name__== "__main__":
    
    parser = argparse.ArgumentParser(description='Input/Output files and other parameters')

    parser.add_argument("--input", dest='input',
                        help ='h5ad object where some QC is going to be applied')
    
    parser.add_argument("--bt_thres", dest='bt_thres', type = float,
                        help ='minimum proportion of cells required for a BATCH to be present')
        
    parser.add_argument("--ct_thres", dest='ct_thres', type = float,
                        help ='minimum proportion of cells required for a CELL TYPE to be present')
    
    parser.add_argument("--min_genes", dest='min_genes', type = int,
                        help ='minimum number of genes to be expressed in a cell')
        
    parser.add_argument("--min_cells", dest='min_cells', type = int,
                        help ='minimum number of cells for a gene to be expressed in')
    
    parser.add_argument('--output', dest='output',
                        help='QCded h5ad object')

    args = parser.parse_args()

    read_h5ad(args.input)
