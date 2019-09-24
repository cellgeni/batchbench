#!/usr/bin/env python3
import numpy as np
import pandas as pd
import scanpy.api as sc
import scanorama
import anndata
import os
import scipy
import sys

import argparse

def save_h5ad(dataset):
    dataset.write(args.output)
    
def correct_scanorama(datasets):
    ''' This function runs Scanorama and saves the corrected object to h5ad file'''
    #Scanorama
    corrected = scanorama.correct_scanpy(datasets)
    #merge Scanorama corrected object
    corrected_dataset = corrected[0].concatenate(corrected[1:], join="inner", batch_key = 'Batch')
    print("Scanorama worked!")
    #save Scanorama corrected object
    save_h5ad(corrected_dataset)


    
def subset_by_batches(dataset):
    #subset the dataset by batches
    N_batches = len(dataset.obs['Batch'].astype('category').cat.categories)
    datasets = []
    for i in range(0, N_batches):
        batch = dataset[dataset.obs['Batch'] == dataset.obs['Batch'].cat.categories[i]]
        datasets.append(batch)

    correct_scanorama(datasets)    

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
