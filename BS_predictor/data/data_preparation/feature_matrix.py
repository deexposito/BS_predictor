''' 
This script is used to create a set of features from a list of 300 binding site structures extracted from the BioLip database.
The resulting feature matrix includes information from 300 pdb structures 
'''

import os
import pandas as pd
from data.classes import Protein, Solvent, Properties, Interactions

# Set the directory to the training PDB files
train_dir = 'data/train_pdbs'

matrices = []

for file in os.listdir(train_dir):
    file_path = os.path.join(train_dir, file)

    if os.path.isfile(file_path):
        prot = Protein(file_path)
        try:
            Solvent().compute_acc(file_path, prot)
            Properties('data').compute_properties(prot)
            Interactions('data').compute_interactions(prot)
        except ValueError:
            continue

        matrices.append(prot.matrix)

full_matrix = pd.concat(matrices)
full_matrix.to_csv('full_matrix.csv', index=False)
