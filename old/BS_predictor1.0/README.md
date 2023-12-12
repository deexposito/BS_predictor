README - FILE

BS_predictor: a machine learning strategy to predict binding sites of an input pdb.

------------------------------------
Denis Expósito Navarro, Marta Alonso Caubilla, Yolanda Andrés López
    Universidad Pompeu Fabra,
    Barcelona, Spain

VERSION 1.0.0  April 2023.

BS_predictor is a supervised learning method for the prediction of binding site 
residues based on the properties and interactions of residues of a protein structure. 
he program is written in python and should work on UNIX like environment.

usage
-----
To work with this model, it is necessary to download all the files provided: 
setup.py
__init__.py
BS_predictor.py
data directory, which contains:
- classes.py
- atom_types.csv
- train_pdbs directory, which contains 300 PDB files


Moreover, note that it is also necessary to install DSSP (Define Secondary 
Structure of Proteins) for the program to work properly. In Linux this can 
be done using the following commands on the terminal:

    sudo apt-get install dssp
    sudo ln -s /usr/bin/mkdssp /usr/bin/dssp

There are two possible ways to run the program:

1. Installing the program BS_predictor: (?)
Run setup.py in the command line:

    python setup.py install

This installs the program BS_predictor.py, making the command BS_predictor
available in the command line.

Run BS_predictor in the command line using the following syntax:

    BS_predictor -p <input_PDB_file> -o <output_PDB_file_name>


2. Running the BS_predictor.py Python script: 

The BS_predictor.py Python script is used to run the program. Simply type 
“python3 BS_predictor.py” to launch the program, and then specify the options 
-p followed by the name of a valid PDB file, and -o followed by the name of 
an output file to save the result, as shown below:

    python3 BS_predictor.py -p <input_PDB_file> -o <output_PDB_file_name>


When using this way of running BS_predictor, it is also necessary to have 
the following packages and versions downloaded:
- pandas >= 2.0.0
- biopython >= 1.81
- networkx >= 3.1
- scipy >= 1.10.1
- sklearn (scikit-learn) >= 1.2.2
- numpy >= 1.24.2

OUTPUT
--------------------
The output is:
- List of residues that belong to the predicted binding site.
- PDB file with information for residues in the predicted binding site to be represented in Chimera or Pymol.

EXAMPLE
--------------------

For the protein structure 4ins.pdb the following residues in the predicted binding site
are obtained:

4ins_A_GLY_1
4ins_A_GLN_5 
4ins_B_GLY_23 
4ins_C_GLY_1 
4ins_D_GLY_23 
