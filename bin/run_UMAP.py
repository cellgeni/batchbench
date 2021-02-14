#!/usr/bin/env python3

#Compute UMAP dimensionality reduction over AnnData objects

import numpy as np
import pandas as pd
import scanpy as sc
import argparse

def save_umap(umap_df, out_name):
    umap_df.to_csv(out_name, header = True, index = True)

def umap_counts_mat(dataset):
    """Compute UMAP over the expression matrix of batch corrected objects, 
    This includes: mnnCorrect, limma, ComBat, Seurat_v3 and Scanorama methods."""
    sc.tl.pca(dataset, n_comps=args.n_pcs)
    sc.pp.neighbors(dataset, n_neighbors=args.n_neighbours, use_rep = 'X_pca')
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"], index = dataset.obs_names,  columns = ["UMAP_1", "UMAP_2"])
    return umap_df

def umap_embedding(dataset):
    """Compute UMAP over the corrected embedding of batch corrected objects,
    This includes: fastMNN and Harmony methods."""
    sc.pp.neighbors(dataset, n_neighbors=args.n_neighbours, use_rep = 'X_corrected_embedding')
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"], index = dataset.obs_names, columns = ["UMAP_1", "UMAP_2"])
    return umap_df

def umap_bbknn(dataset):
    """Compute UMAP over the graph of batch corrected object,
    This includes: BBKNN method."""
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"], index = dataset.obs_names, columns = ["UMAP_1", "UMAP_2"])
    return umap_df

def distribute_datasets(dataset):
    try:
        neighbors = dataset.uns['neighbors']
        print('BBKNN corrected object!')
        save_umap(umap_bbknn(dataset), out_name = args.output_umap)

    except KeyError:
        
        try: 
            dataset.obsm['X_corrected_embedding']
            print('PCA space batch corrected object!')
            save_umap(umap_embedding(dataset), out_name = args.output_umap)

        except KeyError:
            print('Counts matrix batch corrected object!')
            save_umap(umap_counts_mat(dataset), out_name = args.output_umap)


def read_h5ad(args):
	dataset = sc.read(args.input_object)
	distribute_datasets(dataset)
# args
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input_object", 
			dest='input_object',
			type=str,
                        help ='Input h5ad object')
    parser.add_argument("--n_neighbours", 
			dest='n_neighbours', 
			type= int,
			default = 30,
                        help ='Number of nearest neighbours per batch to perform graph correction.')
    parser.add_argument("--n_pcs", 
			dest='n_pcs', 
			type = int,
			default= 25, 
                        help ='Number of PCs for PCA prior to graph correction')
    parser.add_argument("--dim_umap", 
			dest='dim_umap', 
			type = int,
			default= 2, 
                        help ='Number of UMAP components to output')
    parser.add_argument('--output_umap', 
			dest='output_umap',
			type=str, 
                        help='Csv with UMAP coordinate values')
    args = parser.parse_args()

    read_h5ad(args)
