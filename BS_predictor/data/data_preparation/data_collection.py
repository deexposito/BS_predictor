# Import required modules
import pandas as pd
from Bio.PDB import PDBList
import os, shutil 

# Load non-redundant data from the BioLip database
url = "https://zhanggroup.org/BioLiP/download/BioLiP_nr.txt.gz"
df = pd.read_csv(url, sep='\t', header=None, compression='gzip')

# Remove duplicates and extract the first 450 protein interactions to extract information for the model (you can modify the code other interactions)
df = df.drop_duplicates(subset=0)
df_450 = df.head(450)

# Keep only the relevant columns and rename them
df_model = df_450[[0, 1, 4, 5, 8, 9, 19, 20]]
df_model.columns = ['pdb_id', 'receptor_chain', 'ligand_id', 'ligand_chain', 'binding_site_residues',
                    'catalytic_site_residues', 'ligand_residues_seq_num', 'receptor_sequence']

# Selected information for the model
relevant_cols = ['pdb_id', 'receptor_chain', 'ligand_id', 'binding_site_residues']
biolip_relevant = df_model.loc[:, relevant_cols]

# Create a dictionary with information of the binding sites of our selected protein interactions
biolip_dict = {}
for _, row in biolip_relevant.iterrows():
    pdb_id = row["pdb_id"]
    receptor_chain = row["receptor_chain"]
    ligand_id = row["ligand_id"]
    binding_site_residues = row["binding_site_residues"]
    
    if pdb_id not in biolip_dict:
        biolip_dict[pdb_id] = []
        
    biolip_dict[pdb_id].append((receptor_chain,ligand_id,binding_site_residues))

### DOWNLOAD PDB FILES FROM BIOLIP ####
# Download the entire PDB files of the BioLip database
pdbl = PDBList()
for pdb_id in biolip_dict.keys(): # Take each key of the biolip_dict (PDB code)
    pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir="./ent_pdbs")

# Extract binding site information from a dictionary of protein structures and sequences, and store the information in a list for further processing
aa_letter_codes = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
            'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K':
            'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
            'Y': 'TYR', 'V': 'VAL'}

result_bindingsite_list = []
for key, values in biolip_dict.items():
    for value in values:
        chain, _, sequence = value # Store the chain and sequence in a tuple
        
        for aa_code in sequence.split():
            if aa_code[0] in aa_letter_codes.keys():
                aa_name = aa_letter_codes[aa_code[0]]
                residue_num = aa_code[1:]

                result_bindingsite_list.append(f'pdb{key}_{chain}_{aa_name}_{residue_num}')

# Convert the list to a DataFrame
df = pd.DataFrame(result_bindingsite_list, columns=['BSite'])

# Save the DataFrame to a CSV file
df.to_csv('result_bindingsite_list.csv', index=False)


# PDB files were downloaded in .ent format by default, so we transform them into .pdb again
ent_folder = './ent_pdbs' # Folder of .ent files
pdb_folder = './train_pdbs' # Folder where .pdb files will be saved

# Create the folder if it does not exist
if not os.path.exists(pdb_folder):
    os.mkdir(pdb_folder)
    
# Iterate on the .ent files
for file in os.listdir(ent_folder):
    if file.endswith('.ent'):
        ent_path = os.path.join(ent_folder, file)
        pdb_path = os.path.join(pdb_folder, file.replace('.ent', '.pdb'))
        
        # Convert .ent to .pdb
        with open(ent_path, 'r') as f_ent, open(pdb_path, 'w') as f_pdb:
            shutil.copyfileobj(f_ent, f_pdb)