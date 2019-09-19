#!/usr/bin/env python3

#modules
import numpy as np
import pandas as pd
import scanpy.api as sc
import os
import sys

import argparse

def save_umap(umap_df, out_name):
    umap_df.to_csv(out_name, header = True, index = False)

def umap_scanorama(dataset):
    """compute umap coordinates over scanorama object"""
    sc.pp.neighbors(dataset, n_neighbors=30)
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"])
    return umap_df

def umap_bbknnn_after(dataset):
    """compute umap coordinates over bbknn object"""
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"])
    return umap_df

def umap_bbknnn_before(dataset):
    """compute umap coordinates over uncorrected counts of object"""
    sc.pp.neighbors(dataset, n_neighbors=30)
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"])
    return umap_df

def distribute_datasets(dataset):
    try:
        neighbors = dataset.uns['neighbors']
        print('BBKNN corrected object!')

        save_umap(umap_bbknnn_after(dataset), out_name = args.output_bbknn)

    except KeyError:

        print('Scanorama corrected object!')
        save_umap(umap_scanorama(dataset), out_name = args.output_scano)

#read file
def read_h5ad(args):
    #read input h5ad
    dataset = sc.read(args.input)
    print("File read!")

    distribute_datasets(dataset)

#this is the main entry point to the compiler to go through when reading the script
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='Batch corrected object --> .h5ad file')

    parser.add_argument('--output_bbknn', dest='output_bbknn',
                        help='UMAP coordinates of a BBKNN batch corrected object --> .csv file')
    parser.add_argument('--output_scano', dest='output_scano',
                        help='UMAP coordinates of a batch correctec object by other method --> .csv file')
    args = parser.parse_args()

    read_h5ad(args)
