#!/usr/bin/env python3
import numpy as np
import pandas as pd
import scanpy as sc
import scanorama
import anndata
import os
import scipy
import sys

import argparse

def save_h5ad(dataset):
    dataset.write(args.output)
    
def correct_scanorama(dataset_list, cell_metadata):
    ''' This function runs Scanorama and saves the corrected object to h5ad file'''
    #Scanorama
    corrected = scanorama.correct_scanpy(dataset_list)
    #merge Scanorama corrected object
    corrected_dataset = corrected[0].concatenate(corrected[1:], join="inner", batch_key = 'Batch')
    print("Scanorama worked!")
    #append metadata
    corrected_dataset.obs = cell_metadata

    save_h5ad(corrected_dataset)


    
def subset_by_batches(dataset):
    #subset the dataset by batches
    batch_names = dataset.obs['Batch'].astype('category').cat.categories
    N_batches = len(batch_names)
    #save metadata
    cell_metadata = dataset.obs
    #subset dataset by batches
    dataset_list = [dataset[dataset.obs['Batch'] == batch_names[i]] for i in range(0, N_batches)]
    correct_scanorama(dataset_list, cell_metadata)    

def read_h5ad(file):
    dataset = sc.read(file)
    
    subset_by_batches(dataset)
    
#this is the main entry point to the compiler to go through when reading the script
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument("--input", dest='input',
                        help ='h5ad object over which Scanorama is going to be ran')

    parser.add_argument('--output', dest='output',
                        help='Scanorama corrected object.Scanorama corrects the expression matrix')
    args = parser.parse_args()	
    read_h5ad(args.input)
