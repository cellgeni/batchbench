#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scanorama
import anndata
import scipy


def save_h5ad(dataset):
	dataset.write(args.output)
	print("Scanorama done!")
# Scanorama
def correct_scanorama(dataset_list, cell_metadata):
	''' This function runs Scanorama and saves the corrected object to h5ad file'''
	# run Scanorama
	corrected = scanorama.correct_scanpy(dataset_list)
	# merge Scanorama corrected object
	corrected_dataset = corrected[0].concatenate(corrected[1:], 
							join='inner', 
							batch_key = args.batch_key)
	# append metadata
	corrected_dataset.obs = cell_metadata
	
	save_h5ad(corrected_dataset)


# Subset dataset by batches    
def subset_by_batches(dataset):
	batch_names = dataset.obs[args.batch_key].astype('category').cat.categories
	N_batches = len(batch_names)
	# metadata
	cell_metadata = dataset.obs
	# subset dataset by batches
	dataset_list = [dataset[dataset.obs[args.batch_key] == batch_names[i]] for i in range(0, N_batches)]
    
	correct_scanorama(dataset_list, cell_metadata)    

def read_h5ad(file):
	dataset = sc.read(file)
	
	subset_by_batches(dataset)
    
# args
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='h5ad object over which Scanorama is going to be ran')
    parser.add_argument("--batch_key", dest='batch_key',
                        help ='Cell key defining Batch')
    parser.add_argument('--output', dest='output',
                        help='Scanorama corrected object.Scanorama corrects the expression matrix')
    args = parser.parse_args()	
    read_h5ad(args.input)
