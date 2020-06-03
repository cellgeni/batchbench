#!/usr/bin/env python3

#Compute entropy over AnnData objects

import argparse
import numpy as np
import pandas as pd
import scanpy as sc 
import anndata
import scipy
from math import log


def shannon_entropy (x, b_vec, N_b):
    
    tabled_values = b_vec[x > 0].value_counts()/ len(b_vec[x >0]) #class 'pandas.core.series.Series'

    tabled_val = tabled_values.tolist() 
    
    entropy = 0.0
    for element in tabled_val:
        if element != 0:
            entropy += element * log(element)
            
    entropy /= log(N_b)

    return(-entropy) #the entropy formula is the -sum, this is why we include the minus sign

def save_file_to_csv(results):
    results.to_csv(args.output_entropy, header = True, index = False)

def compute_entropy(df, **kwargs):
    #apply function
    batch_entropy = df.apply(shannon_entropy, axis=0, args=(kwargs['batch_vector'],kwargs['N_batches']))
    cell_type_entropy = df.apply(shannon_entropy, axis=0, args=(kwargs['cell_type_vector'] ,kwargs['N_cell_types']))
    print("Entropy calculated!")
    
    results = {'Batch_entropy': batch_entropy, "Cell_type_entropy":cell_type_entropy}
    results = pd.concat(results, axis = 1, keys = ['Batch_entropy', 'Cell_type_entropy'])

    save_file_to_csv(results)

def distribute_datasets(dataset):
    kwargs = {}
    #batch vector(batch id of each cell)
    kwargs['batch_vector'] = dataset.obs[args.batch_key]
    #modify index of batch vector so it coincides with matrix's index
    kwargs['batch_vector'].index = range(0,len(kwargs['batch_vector']))
    #number of batches
    kwargs['N_batches'] = len(dataset.obs['Batch'].astype('category').cat.categories)

    #cell_type vector( betch id of each cell)
    kwargs['cell_type_vector'] = dataset.obs[args.celltype_key]
    #modify index of cell_type vector so it coincides with matrix's index
    kwargs['cell_type_vector'].index = range(0,len(kwargs['cell_type_vector']))
    #number of cell_types
    kwargs['N_cell_types'] = len(dataset.obs['cell_type1'].astype('category').cat.categories)    
    

    try:
        knn_graph = dataset.uns['neighbors']
        print('BBKNN corrected object!') 
    
    except KeyError:
        #Both: pre corrected logcounts and Scanorama counts enter through this way.
        #compute neighbors
        sc.tl.pca(dataset, n_comps = args.n_pcs)
        sc.pp.neighbors(dataset, n_neighbors = args.n_neighbors, knn = True)
                        
    #knn graph
    knn_graph = dataset.uns['neighbors']['connectivities']
    #transforming csr_matrix to dataframe
    df = pd.DataFrame(knn_graph.toarray())
    
    compute_entropy(df, **kwargs)

def read_h5ad(dataset):
    
    #read input h5ad
    return sc.read(dataset)
    print("File read!")
    

# args
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input_object", 
			dest='input_object',
			type=str,
                        help ='Input h5ad object')
    parser.add_argument("--batch_key", 
			dest='batch_key',
                        type=str,
			default='Batch',
			help ='Cell key defining Batch')
    parser.add_argument("--celltype_key", 
			dest='celltype_key', 
			type = str,
			default= 'cell_type1', 
			help ='Cell key defining cell type')
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
    parser.add_argument('--output_entropy', 
			dest='output_entropy',
			type=str, 
                        help='Csv with entropy values')
    args = parser.parse_args()

    read_h5ad(args)
