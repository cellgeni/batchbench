#modules
import numpy as np
import pandas as pd
import scanpy.api as sc #should it be .API?

import anndata
import os
import scipy
import sys

from math import log

import argparse



def shannon_entropy (x, b_vec, N_b):
    
    tabled_values = b_vec[x > 0].value_counts()/ len(b_vec[x >0]) #class 'pandas.core.series.Series'

    tabled_val = tabled_values.tolist() 
    
    entropy = 0.0
    for element in tabled_val:
        if element != 0:
            entropy += element * log(element)
            
    entropy /= log(N_b)

    return(-entropy) #the entropy formula is the -sum, this is why we include the minus sign

def compute_entropy(dataset):
    
    try:
        knn_graph = dataset.uns['neighbors']
        
    except KeyError:
        #compute neighbors
        sc.pp.neighbors(dataset, n_neighbors = 30, knn = True)
                        
    #knn graph
    knn_graph = dataset.uns['neighbors']['connectivities']
    #transforming csr_matrix to dataframe
    df = pd.DataFrame(knn_graph.toarray())
    
    #batch vector( betch id of each cell)
    batch_vector = dataset.obs['batch']
    #modify index of batch vector so it coincides with matrix's index
    batch_vector.index = range(0,len(batch_vector))
    #number of batches
    N_batches = len(dataset.obs['batch'].astype('category').cat.categories)
    
    #cell_type vector( betch id of each cell)
    cell_type_vector = dataset.obs['cell_type']
    #modify index of cell_type vector so it coincides with matrix's index
    cell_type_vector.index = range(0,len(cell_type_vector))
    #number of cell_types
    N_cell_types = len(dataset.obs['cell_type'].astype('category').cat.categories)
    
    #apply function
    global_entropy = df.apply(shannon_entropy, axis=0, args=(batch_vector, N_batches))
    cell_type_entropy = df.apply(shannon_entropy, axis=0, args=(cell_type_vector, N_cell_types))
    print("Entropy calculated!")
    
    results = {'global_entropy':global_entropy, "cell_type_entropy":cell_type_entropy}
    results = pd.concat(results, axis = 1).reset_index()
    
    #save pandas df to output specified path
    results.to_csv(args.output, header = True)
    

def read_h5ad(args):
    
    #read input h5ad
    dataset = sc.read(args.input)
    print("File read!")
    
    compute_entropy(dataset)

#this is the main entry point to the compiler to go through when reading the script

if __name__== "__main__":
    
    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument('input',
                        help='batch corrected file')
    
    parser.add_argument('output',
                        help='object with entropy value per cell')

    args = parser.parse_args()
    
    read_h5ad(args)