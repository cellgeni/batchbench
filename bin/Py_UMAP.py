#!/usr/bin/env python3
import numpy as np
import pandas as pd
import scanpy as sc
import argparse

def save_umap(umap_df, out_name):
    umap_df.to_csv(out_name, header = True, index = False)

def umap_counts_mat(dataset):
    """Compute UMAP over the expression matrix of batch corrected objects, 
    This includes: mnnCorrect, limma, ComBat, Seurat_v3 and Scanorama methods."""
    sc.tl.pca(dataset, n_comps= 20)
    sc.pp.neighbors(dataset, n_neighbors=30, use_rep = 'X_pca')
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"], columns = ["UMAP_1", "UMAP_2"])
    return umap_df

def umap_embedding(dataset):
    """Compute UMAP over the corrected embedding of batch corrected objects,
    This includes: fastMNN and Harmony methods."""
    sc.pp.neighbors(dataset, n_neighbors=30, use_rep = 'X_corrected_embedding')
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"], columns = ["UMAP_1", "UMAP_2"])
    return umap_df

def umap_bbknn(dataset):
    """Compute UMAP over the graph of batch corrected object,
    This includes: BBKNN method."""
    sc.tl.umap(dataset, n_components=2)
    umap_df = pd.DataFrame(dataset.obsm["X_umap"], columns = ["UMAP_1", "UMAP_2"])
    return umap_df

#def umap_bbknnn_before(dataset):
#    """compute umap coordinates over uncorrected counts of object"""
#    sc.pp.neighbors(dataset, n_neighbors=30)
#    sc.tl.umap(dataset, n_components=2)
#    umap_df = pd.DataFrame(dataset.obsm["X_umap"])
#    return umap_df

def distribute_datasets(dataset):
    try:
        neighbors = dataset.uns['neighbors']
        print('BBKNN corrected object!')
        save_umap(umap_bbknn(dataset), out_name = args.output)

    except KeyError:
        
        try: 
            dataset.obsm['X_corrected_embedding']
            print('PCA space batch corrected object!')
            save_umap(umap_embedding(dataset), out_name = args.output)

        except KeyError:
            print('Counts matrix batch corrected object!')
            save_umap(umap_counts_mat(dataset), out_name = args.output)

#def save_corrected_h5ad(dataset):
#	dataset.write(args.output_h5ad)

def read_h5ad(args):
	dataset = sc.read(args.input)
	distribute_datasets(dataset)
#	save_corrected_h5ad(dataset)

if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='Batch corrected object --> .h5ad file')

    parser.add_argument('--output', dest='output',
                        help='UMAP coordinates of a batch corrected object --> .csv file')
    
 #   parser.add_argument('--output_h5ad', dest='output_h5ad',
  #                      help='h5ad batch corrected object, saved for additional analysis --> .h5ad file')
    
    args = parser.parse_args()

    read_h5ad(args)
