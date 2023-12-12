#!/usr/bin/env python
# coding: utf-8
### Authors: Marta Alonso, Yolanda Andres and Denis Exposito ###

# Required imports
from Bio import PDB
from Bio.PDB import PDBParser, PDBList, is_aa, NeighborSearch, DSSP, PDBIO
from Bio.SeqUtils import IUPACData
import os, sys, getopt, warnings, shutil, scipy
import pandas as pd
from data.classes import Input_Protein, Protein, Neighborhood, Solvent, Properties, Interactions, Prediction, result_bindingsite_list, get_residue_index, residues_to_keep, read_atom_types
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.utils import resample
import numpy as np
import networkx as nx

# Obtain current directory and create a dictionary with the directory for the directories data and train_pdbs
base_dir = os.path.dirname(__file__)
config = {'data_dir': os.path.join(base_dir,'data/'), 'train_dir': os.path.join(base_dir,'data/train_pdbs')}

# Main function of the binding site predictor program
def main(argv):
    
    '''Take a list of command line arguments (paths to input protein structure file and output directory) as input and predict the binding sites of the input protein, returning a list of residues and a new PDB file containing them.'''
    
    # Parse the arguments and obtain the paths to the input protein structure file and the output directory
    try:
        opts, args = getopt.getopt(argv, "p:o:", ['pdb_file=','out='])
    except getopt.GetoptError: # If the input is incorrect, print an error message
        print ('Incorrect input. Please use the following syntax: BS_predictor.py -p <pdb_file> -o <out_dir>')
        sys.exit(2)

    protein_file = ''
    out_dir = ''

    for opt, arg in opts:
        if opt in ("-p", "--pdb_file"):
            protein_file = arg
        elif opt in ("-o", "--out"):
            out_dir = arg
    
    # Instantiate an Input_Protein object with the path to the input protein structure file
    target_protein = Input_Protein(protein_file)
    
    # Compute the solvent accessibility of the input protein 
    print(f"Computing exposure information of protein {target_protein.pdb_id}...")
    Solvent().compute_acc(protein_file, target_protein)
    
    # Compute the residue properties of the input protein
    print(f"Computing residue properties of protein {target_protein.pdb_id}...")
    Properties(config['data_dir']).compute_properties(target_protein)
    
    # Compute the protein interactions of the input protein
    print(f"Computing interactions of protein {target_protein.pdb_id}...")
    Interactions(config['data_dir']).compute_interactions(target_protein)
    
    print(f"Predicting the binding sites of protein {target_protein.pdb_id}...")
    
    # Generate a matrix of features (full_matrix) for the training PDB files in the directory data/train_pdbs
    matrices = []
    
    for file in os.listdir(config['train_dir']):
        file_path = os.path.join(config['train_dir'], file)
    
        if os.path.isfile(file_path):
            prot = Protein(file_path)
            try:
                Solvent().compute_acc(file_path, prot)
                Properties(config['data_dir']).compute_properties(prot)
                Interactions(config['data_dir']).compute_interactions(prot)
            except ValueError:
                continue
        
            matrices.append(prot.matrix)

    full_matrix = pd.concat(matrices)
    
    # Use this feature matrix (full_matrix) to train a Prediction object and use it to predict the binding sites of the input protein
    BS_prediction = Prediction(full_matrix, target_protein, out_dir).predict_binding_sites()
    
    # Generate and print the list of predicted binding sites and create an output PDB file containing these residues
    print(f"The list of predicted binding sites of protein {target_protein.pdb_id} is the following:")
    out_list = Prediction(full_matrix, target_protein, out_dir).get_output(BS_prediction)
    print(out_list)
    
    print(f"A new pdb file has been generated with the predicted binding sites of protein {target_protein.pdb_id}.")

    print("Program finished correctly.")

# Check if the script is being run as the main program and, if so, call the main() function passing the command line arguments
if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Temporarily suppress warnings while running the main() function
        main(sys.argv[1:])
