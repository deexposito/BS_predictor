#!/usr/bin/env python
# coding: utf-8
### Authors: Marta Alonso, Yolanda Andres and Denis Exposito ###

# Required imports
from Bio.PDB import PDBParser, PDBList, is_aa, NeighborSearch, DSSP, PDBIO
from Bio.SeqUtils import IUPACData
import os, sys, getopt, warnings, shutil, scipy, joblib
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.utils import resample
import numpy as np
import networkx as nx
from .data.binding_sites import result_bindingsite_list

# Obtain current directory and create a dictionary with the directory for the directories data and train_pdbs
base_dir = os.path.dirname(__file__)
config = {'data_dir': os.path.join(base_dir,'data/'), 'model': os.path.join(base_dir, 'model/BSmodel.pkl')}

# List of amino acids 
residues_to_keep = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

#### PROTEIN CLASS #### (Used for generating a feature matrix of training data)
class Protein:
    '''Receive the path of a PDB file as a string and return a "Protein" object containing information of the protein.'''

    def __init__(self, file_name):
        '''Initialize the protein matrix object that contains information about the residues and whether each of them belong to a binding site (1) or not (0). Does not return any output.'''

        self.file_name = file_name
        self.pdb_id = os.path.basename(self.file_name).split('.')[0]
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure(self.pdb_id, self.file_name)
        self.init_matrix()

        # Neighborhood radios = 5A
        self.neighborhood = Neighborhood().run(self, 5)
        
    def init_matrix(self):
        '''Initialize the protein matrix object that contains information about the residues. Does not return any output.'''

        # List all residues
        residue_list = list(self.structure[0].get_residues())
        residue_names = []
        for residue in residue_list:
            if (is_aa(residue) or residue in residues_to_keep) and residue.id[0] == ' ': 
                # Residue id: pdbID_chain_resName_resNumber
                residue_name = get_residue_index(residue)
                residue_names.append(residue_name)

        # Create empty matrix
        self.matrix = pd.DataFrame(index=residue_names)
        self.matrix.index.name = 'res_name'
        
        # add a new column "binding_site" with predetermined value 0
        self.matrix['binding_site'] = 0

        # Use the function isin() to determine if a residue is in the list of binding sites obtained from the BioLip database
        is_in_list = self.matrix.index.isin(result_bindingsite_list)

        # Assign value of 1 to the column "binding_site" to the residues that belong to a binding site (according to BioLip)
        self.matrix.loc[is_in_list, 'binding_site'] = 1

# Define a function to obtain the index of a residue 
def get_residue_index(residue):
    '''Generate a unique index for a given residue in a PDB structure, based on its ID, chain, name and number.'''

    pdb_id = residue.get_full_id()[0]
    chain = residue.get_parent().get_id()
    name = residue.get_resname()
    num = str(residue.get_id()[1])

    return pdb_id + '_' + chain + '_' + name + '_' + num


#### PROTEIN CLASS FOR INPUT FILE ####
# Define the class Input_Protein to handle protein structure information of the input protein
class Input_Protein:
    '''Receive the path of an input PDB file as a string and return a "Input_Protein" object containing information of the protein.'''
    
    # Initialize the attributes of the object
    def __init__(self, file_name):
        self.file_name = file_name
        self.pdb_id = os.path.basename(self.file_name).split('.')[0]
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure(self.pdb_id, self.file_name)
        self.init_matrix()
        self.neighborhood = Neighborhood().run(self, 5) # Neighborhood radios = 5A
     
    # Initialize the protein matrix object
    def init_matrix(self):
        '''Initialize the protein matrix object that contains information about the residues. Does not return any output.'''

        # List all residues
        residue_list = list(self.structure[0].get_residues())
        residue_names = []
        for residue in residue_list:
            # Check if the residue is a valid amino acid
            if (is_aa(residue) or residue in residues_to_keep) and residue.id[0] == ' ':
                # Check if residue has all heavy atoms (except H) and, if not, raise a ValueError
                try:
                    heavy_atoms = []
                    for atom_name in ['N', 'CA', 'C', 'O']:
                        if residue.has_id(atom_name):
                            heavy_atoms.append(residue[atom_name])
                        else:
                            raise ValueError("PDB file only contains information about C-alpha atoms.")
                    if not all(atom.name in [a.name for a in residue.get_unpacked_list() if a.element != 'H'] for atom in heavy_atoms):
                        raise ValueError("PDB file only contains information about C-alpha atoms.")
                except ValueError:
                    raise
                    
                # Residue id: pdbID_chain_resName_resNumber
                residue_name = get_residue_index(residue)
                residue_names.append(residue_name)

        # Initialize empty matrix with the residue names as indices
        self.matrix = pd.DataFrame(index=residue_names)
        self.matrix.index.name = 'res_name'


#### NEIGHBORHOOD CLASS ####

# Define class Neighborhood to get the neighbor residues of all atoms
class Neighborhood:
    '''Compute the neighbor residues for a given protein and distance threshold. Store neighbor residues in a dictionary.'''

    def run(self, protein, distance):
        '''Take a Protein object and a distance threshold and return a dictionary where each residue name is the key and the corresponding value is a set of neighbor residue names.'''

        self.adjacent_list = {}
        # list all atoms to compute neighbor residues
        self.atomlist = list(protein.structure.get_atoms())
        # list all residues
        self.residue_list = list(protein.structure.get_residues())
        # object to compute neighbor residues
        self.n = NeighborSearch(self.atomlist)

        for residue in self.residue_list:
            # Check if the residue is a valid amino acid
            if is_aa(residue):
                # Generate a unique residue index: pdbID_resChain_resName_resNumber
                residues_name = get_residue_index(residue)
                # list of neighbor residues
                neighbor_list = []
                for atom in residue.get_atoms():
                    # Compute neighbors within the specified distance threshold and add them to the list
                    neighbor_list += self.n.search(atom.get_coord(), distance, level='R')

                # Save unique neighbor residue names and delete copies
                neighbor_set = set()
                for neighbor in neighbor_list:
                    if is_aa(neighbor):
                        # pdbID_residueChain_residueName_residueNumber
                        neighbor_name = get_residue_index(neighbor)
                        # Avoid adding the same residue
                        if neighbor_name != residues_name:
                            neighbor_set.add(neighbor_name)
                            
                # Save the set of neighbor residues as the value for the key residues_name in the adjacent_list dictionary
                self.adjacent_list[residues_name] = neighbor_set
        return self.adjacent_list




# Function needed to compute the residue properties and interactions 
def read_atom_types(data_dir):
    '''Read the atom_types.csv file containing information about the atom types and return information in a dictionary where the keys are the atom names and the values are lists of properties for each atom type.'''

    types = {}
    # Open the CSV file located in the directory passed as argument and read each line
    with open(os.path.join(data_dir,"atom_types.csv")) as in_file:
        for line in in_file:
            record = line.strip().split(',')
            atomName = record[0] + '_' + record[1]# Use first 2 elements of the line (element symbol and atom type) as atom name
            # If the line has < 3 elements, skip it since it corresponds to atoms without properties
            if len(record) < 3:
                continue
            else:
                types[atomName] = record[2:] # Store remaining elements of the line (atom properties) as values of the dictionary
    return types



#### ATOM PROPERTIES CLASS ####
# Define the class Properties to compute the properties of atoms and residues in a protein structure
class Properties:
    '''Define a set of methods to compute properties of atoms and residues in a protein structure based on the file atom_types.csv.'''

    # Initialize the attributes of the object
    def __init__(self, data_dir):
        self.types = read_atom_types(data_dir)

    def compute_properties(self, prot):
        '''Take a Protein object and compute the properties of each atom of each residue, adding new columns to the matrix of the Protein object.'''

        parser = PDBParser()

        # Insert new columns in protein matrix
        columns = ['DON','HPB','ACP','NEG','ARM','POS', 'CYS']
        for c in columns:
            prot.matrix[c] = 0

        for residue_name in prot.matrix.index:
            pdb_id, res_chain, res_name, res_number = residue_name.split('_')
            res_number = int(res_number)
            structure = prot.structure
    
            for model in structure:
                for chain in model:
                    for residue in chain:
                        # Verify if the name and number of each residue match the ones in the matrix index
                        if residue.get_resname() == res_name and residue.get_id()[1] == res_number:
                            for atom in residue:
                                atom_name = residue.get_resname() + '_' + atom.get_name()
                        
                                # # Verify if the atom has properties and, if so, increase the corresponding property count in the matrix
                                if atom_name in self.types:
                                    for c in columns:
                                        if c in self.types[atom_name]:
                                            prot.matrix.loc[residue_name, c] += 1
                                            
                            # If the residue is a Cysteine, set the corresponding CYS column to 1.                
                            if res_name == 'CYS':
                                prot.matrix.loc[residue_name,'CYS'] = 1
        return prot.matrix




#### INTERACTIONS CLASS ####
# Define the class Interactions to compute the interactions of atoms and residues in a protein structure
class Interactions:
    '''Define a set of methods to compute the interactions of atoms and residues in a protein structure based on the file atom_types.csv.'''

    # Initialize the attributes of the object
    def __init__(self, data_dir):
        self.atom_types = read_atom_types(data_dir)

    def compute_interactions(self, prot):
        '''Take a Protein object and compute the number and type of interactions of each residue, adding new columns to the matrix of the Protein object.'''

        # Insert new columns in the protein matrix and set values to 0
        columns = ['aromatic_stacking','disulfide_bridge','hydrogen_bond', 'hydrophobic', 'repulsive', 'salt_bridge']
        for c in columns:
            prot.matrix[c] = 0

        # Create a Residue Graph to store information about the interactions between residues
        prot.residue_graph = nx.Graph()

        # Create a dataframe to keep track of which residues have been visited
        residue_info_df = pd.DataFrame(columns = ['visited'], index = prot.matrix.index)
        residue_info_df['visited'] = False

        for residue_name in prot.matrix.index:
            pdb_id, res_chain, res_name, res_number = residue_name.split('_')
            res_number = int(res_number)

            # Verify if the residue exists in the protein data structure
            if res_number in prot.structure[0][res_chain]:
                residue = prot.structure[0][res_chain][res_number]

                # Computes its interactions with its neighboring residues
                for adj_res_name in prot.neighborhood[residue_name]:

                    # Verify if the residue exists in the protein matrix
                    if adj_res_name in prot.matrix.index:

                        # Verify if the residues was visitated
                        if not residue_info_df.loc[adj_res_name, 'visited'].any():

                            # split residue name
                            pdb_id, adj_chain, adj_name , adj_number = adj_res_name.split('_')
                            adj_number = int(adj_number)

                            # Verify if the neighbor residue exixts on protein data structure
                            if adj_number in prot.structure[0][adj_chain]:
                                adj_residue = prot.structure[0][adj_chain][adj_number]

                                # Access residue atoms
                                for atom_1 in residue:
                                    atom1_name = res_name + '_' + atom_1.get_name().strip()
                                    if atom1_name in self.atom_types:

                                        # Access neighbor residue atoms
                                        for atom_2 in adj_residue:
                                            atom2_name = adj_name + '_' + atom_2.get_name().strip()
                                            if atom2_name in self.atom_types:

                                                # Compute the distance between atoms
                                                distance = atom_1 - atom_2
                                        
                                                # Test the different interaction types
                                                # Aromatic stacking: 2 aromatic residues and 1.5-3.5 A
                                                if ('ARM' in self.atom_types[atom1_name]) and ('ARM' in self.atom_types[atom2_name]) and (distance >= 1.5) and (distance <= 3.5):
                                                    prot.matrix.loc[residue_name,'aromatic_stacking'] += 1
                                                    prot.matrix.loc[adj_res_name,'aromatic_stacking'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                
                                                # Hydrogen bond: a donor and an acceptor residue and 2-3 A
                                                if ('ACP' in self.atom_types[atom1_name]) and ('DON' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.0):
                                                    prot.matrix.loc[residue_name,'hydrogen_bond'] += 1
                                                    prot.matrix.loc[adj_res_name,'hydrogen_bond'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('DON' in self.atom_types[atom1_name]) and ('ACP' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.0):
                                                    prot.matrix.loc[residue_name,'hydrogen_bond'] += 1
                                                    prot.matrix.loc[adj_res_name,'hydrogen_bond'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                    
                                                # Hydrophobic interaction: 2 hydrophobic residues and 2-3.8 A
                                                if ('HPB' in self.atom_types[atom1_name]) and ('HPB' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 3.8):
                                                    prot.matrix.loc[residue_name,'hydrophobic'] += 1
                                                    prot.matrix.loc[adj_res_name,'hydrophobic'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                    
                                                # Repulsive interaction: 2 positive or 2 negative residues and 2-6 A
                                                if ('POS' in self.atom_types[atom1_name]) and ('POS' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    prot.matrix.loc[residue_name,'repulsive'] += 1
                                                    prot.matrix.loc[adj_res_name,'repulsive'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('NEG' in self.atom_types[atom1_name]) and ('NEG' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    prot.matrix.loc[residue_name,'repulsive'] += 1
                                                    prot.matrix.loc[adj_res_name,'repulsive'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                
                                                # Salt bridge: a positive and a negative residue and 2-6 A
                                                if ('NEG' in self.atom_types[atom1_name]) and ('POS' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    prot.matrix.loc[residue_name,'salt_bridge'] += 1
                                                    prot.matrix.loc[adj_res_name,'salt_bridge'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                if ('POS' in self.atom_types[atom1_name]) and ('NEG' in self.atom_types[atom2_name]) and (distance >= 2.0) and (distance <= 6.0):
                                                    prot.matrix.loc[residue_name,'salt_bridge'] += 1
                                                    prot.matrix.loc[adj_res_name,'salt_bridge'] += 1
                                                    prot.residue_graph.add_edge(residue_name,adj_res_name)
                                                    
                                                # Disulfide bridge: 2 disulfide bridge forming atoms and 3-7.5 A
                                                if ('SSB' in self.atom_types[atom1_name]) and ('SSB' in self.atom_types[atom2_name]):
                                                    # check if residues form a disulfide bridge
                                                    if abs(residue['CA'] - adj_residue['CA']) <= 4.5:
                                                        prot.matrix.loc[residue_name,'disulfide_bridge'] += 1
                                                        prot.matrix.loc[adj_res_name,'disulfide_bridge'] += 1
                                                        prot.residue_graph.add_edge(residue_name, adj_res_name)

        residue_info_df.loc[residue_name,'visited'] = True

        return prot.matrix


#### SOLVENT CLASS ####
# Define the class Solvent to compute the solvent accessibility of each residue in a protein structur
class Solvent:
    '''Compute the solvent accessibility of each residue in a protein structure.'''
        
    def compute_acc(self, pdb_file, prot):
            '''Take a protein structure (PDB) file and a Protein object and compute the solvent accessibility of each residue using the DSSP algorithm, adding a new column (called "ACC") to the matrix of the Protein object.'''

            structure = prot.structure
            model = structure[0]
            dssp = DSSP(model, pdb_file, dssp='mkdssp') # Compute the solvent accessibility of each residue using the DSSP algorithm
            list_acc = []
            
            for key in dssp.keys():
                try:
                    residue = structure[0][key[0]][key[1]]
                except KeyError:
                    continue
                # Check if the residue is a valid amino acid
                if (is_aa(residue) or residue in residues_to_keep):
                    acc = dssp[(key)][3]
                    list_acc.append(acc)
            # Add the computed values as a new column named "ACC" to the matrix of the Protein object    
            prot.matrix["ACC"] = list_acc
            return prot.matrix



#### PREDICTION CLASS ####
# Define the class Prediction to predict the binding site residues of a protein using a Random Forest Classification approach
class Prediction:
    '''Predict the binding site residues of a target protein from a feature matrix, using a Random Forest classification algorithm.'''

    # Initialize the attributes of the object
    def __init__(self, target_protein, out_dir):
        self.target_protein = target_protein # Protein object
        self.out_dir = out_dir if out_dir else 'predicted_binding_residues.pdb' # Output directory and file name for the PDB file
    
    # Define a method to predict the binding site residues of a protein
    def predict_binding_sites(self):
        '''Predict the binding site residues of a target protein using a Random Forest Classification algorithm.'''

        model = joblib.load(config['model'])
        
        # Predict the binding sites for the target sample
        BS_prediction = model.predict(self.target_protein.matrix)
        
        return BS_prediction
    
    # Define a method to obtain the desired output (list of binding site residues and PDB file)    
    def get_output(self, BS_prediction):
        '''Take a binding site prediction object and return a list of residues predicted as binding site residues of a target protein and a PDB file with these binding site residues.'''

        # Obtain the list of residues predicted to be binding sites
        selected_residues = self.target_protein.matrix.index[BS_prediction == 1].tolist()
        structure = self.target_protein.structure

        # Get the chain and residue information for the selected residues
        selected_residues_info = [(residue_id.split('_')[1], residue_id.split('_')[2], int(residue_id.split('_')[3])) for residue_id in selected_residues]

        # Iterate over all atoms and delete those that are not in the selected list
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                for residue in list(chain.get_residues()):
                    resseq = residue.get_id()[1]
                    resname = residue.get_resname()
                    if (chain_id, resname, resseq) not in selected_residues_info:
                        chain.detach_child(residue.get_id())

        # Save the modified structure to a new PDB file
        io = PDBIO()
        io.set_structure(structure)
        io.save(self.out_dir)
        
        # Return the list of predicted binding site residues
        return selected_residues

