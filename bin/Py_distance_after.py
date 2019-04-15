#!/usr/bin/env python3

#modules
import numpy as np
import pandas as pd
import scanpy.api as sc 

import anndata
import os
import scipy
import sys
from sklearn.metrics.pairwise import cosine_similarity
import argparse

def save_file_to_csv(results):
    results.to_csv(args.output, header = True, index = False)

def merge_results_bbknn(bbknn_results_dict, before_results_dict):
    "Merge dictionaries containing distances after bbknn and distances before bbknn"
    dict_bbknn_bef_aft = {**before_results_dict, **bbknn_results_dict}
    #we construct the dataframe this way(instead of DataFrame.from_dict) so that columns can have different lengths!
    results = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dict_bbknn_bef_aft.items() ]))
    
    save_file_to_csv(results)


def calc_cosine_similarity_scanorama(total_list, **kwargs):
    """Calculate the cosine similarity over the subsets generated in subset_by_scanorama"""
    sim_list = []
    for i in range(0, len(total_list)):
        similarities = cosine_similarity(total_list[i].X)
        sim_list.append(similarities)    
    
    column_names = kwargs['batch_names'] + kwargs['cell_type_names']
    print(column_names) ########
    vals_dict = {}

    for i in range(0, len(sim_list)):
        #get half of the values from the symmetrix matrix, without diagonal(k = 1)
        vals = sim_list[i][np.triu_indices(len(sim_list[i]), k = 1)]
        vals = vals.tolist()
        #name is counts
        if i < kwargs['N_batches']:
            vals_dict["COUNTS_dist_BATCH_" + column_names[i]]  = vals
        else:
            vals_dict["COUNTS_dist_CELLTYPE_" + column_names[i]]  = vals

        #we construct the dataframe this way(instead of DataFrame.from_dict) so that columns can have different lengths!
        results_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in vals_dict.items() ]))
        
    return results_df
    #save_file_to_csv(results_df)

def subset_by_scanorama(dataset, **kwargs):
    """Here we subset the datasets By Batch and By cell type"""
    
    #subset the dataset by batches
    batch_datasets = []
    for i in range(0, kwargs['N_batches']):
        batch = dataset[dataset.obs['Batch'] == dataset.obs['Batch'].cat.categories[i]]
        batch_datasets.append(batch)

    #subset the dataset by batches
    cell_datasets = []
    for i in range(0, kwargs['N_cell_types']):
        cell_type = dataset[dataset.obs['cell_type1'] == dataset.obs['cell_type1'].cat.categories[i]]
        cell_datasets.append(cell_type)

    #merged list of subsetted arrays
    total_list = batch_datasets + cell_datasets
    
    return total_list
    
def distances_after_bbknnn(dist_df, **kwargs):
    """This function subsets the distance graph (symmetric) by batch and cell type,
    and attaches the values of the upper triangle to a dictonary"""
    
    bbknn_results_dict = {}
    #subset by batches
    for i in range(0, kwargs['N_batches'] + kwargs['N_cell_types']):
        if i < kwargs['N_batches']:
            dist_df.index = kwargs['batch_vector']
            dist_df.columns = kwargs['batch_vector']
            subset = dist_df.filter(items = [kwargs['batch_names'][i]],axis = 0).filter(items = [kwargs['batch_names'][i]],axis = 1)
            #turn to numpy array again
            subset = subset.values
            #get values from upper triangle without diagonal
            vals = subset[np.triu_indices(len(subset), k = 1)]
            vals = vals.tolist()
            bbknn_results_dict["BBKNN_Batch_dist_" + kwargs['batch_names'][i]] = vals
    
    else:
        dist_df.index = kwargs['cell_type_vector']
        dist_df.columns = kwargs['cell_type_vector']
        for j in range(0,kwargs['N_cell_types']):
            subset = dist_df.filter(items = [kwargs['cell_type_names'][j]],axis = 0).filter(items = [kwargs['cell_type_names'][j]],axis = 1)
            #turn to numpy array again
            subset= subset.values
            #get values from upper triangle without diagonal
            vals = subset[np.triu_indices(len(subset), k = 1)]
            vals = vals.tolist()
            bbknn_results_dict["BBKNN_Cell_type_dist_" + kwargs['cell_type_names'][j]] = vals
    
    return bbknn_results_dict
       
def distances_before_bbknnn(dataset, **kwargs):
    #compute neighbors
    sc.pp.neighbors(dataset, n_neighbors=30, n_pcs=None, knn = True)
    distances = dataset.uns['neighbors']['distances']
    dist_df = pd.DataFrame(distances.toarray())

    before_results_dict = {}
    #subset by batches
    for i in range(0, kwargs['N_batches'] + kwargs['N_cell_types']):
        if i < kwargs['N_batches']:
            dist_df.index = kwargs['batch_vector']
            dist_df.columns = kwargs['batch_vector']
            subset = dist_df.filter(items = [kwargs['batch_names'][i]],axis = 0).filter(items = [kwargs['batch_names'][i]],axis = 1)
            #turn to numpy array again
            subset= subset.values
            #get values from upper triangle without diagonal
            vals = subset[np.triu_indices(len(subset), k = 1)]
            vals = vals.tolist()
            before_results_dict["Bef_Batch_dist_" + kwargs['batch_names'][i]] = vals
    
    else:
        dist_df.index = kwargs['cell_type_vector']
        dist_df.columns = kwargs['cell_type_vector']
        
        for j in range(0,kwargs['N_cell_types']):
            subset = dist_df.filter(items = [kwargs['cell_type_names'][j]],axis = 0).filter(items = [kwargs['cell_type_names'][j]],axis = 1)
            #turn to numpy array again
            subset= subset.values
            #get values from upper triangle without diagonal
            vals = subset[np.triu_indices(len(subset), k = 1)]
            vals = vals.tolist()
            before_results_dict["Bef_Cell_type_dist_" + kwargs['cell_type_names'][j]] = vals
    
    return before_results_dict
       
def distribute_datasets(dataset):
    #creating keyword arguments for other functions (this helps avoid repeating code)
    kwargs = {}
    kwargs['batch_vector'] = dataset.obs['Batch']
    kwargs['batch_names'] = kwargs['batch_vector'] .astype('category').cat.categories.tolist()
    kwargs['N_batches'] = len(kwargs['batch_names'])

    kwargs['cell_type_vector'] = dataset.obs['cell_type1']
    kwargs['cell_type_names'] =  kwargs['cell_type_vector'].astype('category').cat.categories.tolist()
    kwargs['N_cell_types'] = len(kwargs['cell_type_names'])
    
    print(dataset)
    try:
        neighbors = dataset.uns['neighbors']
        distances = neighbors['distances']
        dist_df = pd.DataFrame(distances.toarray())
        print('BBKNN corrected object!')
                                                                                                     
        merge_results_bbknn(distances_after_bbknnn(dist_df, **kwargs), distances_before_bbknnn(dataset, **kwargs))
        
    
    except KeyError:
        
        print('Scanorama corrected object!')
        save_file_to_csv(calc_cosine_similarity_scanorama(subset_by_scanorama(dataset, **kwargs), **kwargs))         
        
def read_h5ad(x):    
    #read input h5ad
    return sc.read(x)
  

if __name__== "__main__":
    
    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='h5ad object over which cosine similarities are going to be calculated')
    
    parser.add_argument('--output', dest='output',
                        help='cosine similarities for batch and cell type subsets stored in a dictionary')

    args = parser.parse_args()

    distribute_datasets(read_h5ad(args.input))
    
