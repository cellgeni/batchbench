#!/usr/bin/env python3
  
#modules
import numpy as np
import pandas as pd
import scanpy.api as sc
import scanorama
import anndata
import os
import scipy
import sys

import argparse

def correct(datasets):
    ''' This function runs Scanorama and saves the corrected object to h5ad file'''
    #Scanorama
    corrected = scanorama.correct_scanpy(datasets)

    #merge Scanorama corrected object
    corrected_dataset = corrected[0].concatenate(corrected[1:], join="inner", batch_key = 'Batch')
    print("Scanorama worked!")

    #save Scanorama corrected object
    corrected_dataset.write(args.output)
    print("Corrected object saved!")



def main(args):
    '''This function reads an h5ad file, splits the file into its batches and calls the correct function '''
    #read input h5ad
    dataset = sc.read(args.input)
    print("File read!")

    #subset the dataset by batches
    N_batches = len(dataset.obs['Batch'].astype('category').cat.categories)
    datasets = []
    for i in range(0, N_batches):
        batch = dataset[dataset.obs['Batch'] == dataset.obs['Batch'].cat.categories[i]]
        datasets.append(batch)

    print(datasets)
    correct(datasets)

#this is the main entry point to the compiler to go through when reading the script
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files')

    parser.add_argument('input',
                        help='input file to be batch corrected')

    parser.add_argument('output',
                            help='output file batch corrected')

    args = parser.parse_args()
	
    main(args)
