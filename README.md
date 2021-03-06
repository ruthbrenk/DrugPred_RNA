DrugPred_RNA
by Illimar Rekand, Ruth Brenk
University of Bergen, 2020


To run DrugPred_RNA, you need Python 3.6 and R version 4.0.3 and UCSF DOCK 3.6 (http://dock.compbio.ucsf.edu/)

The following Python packages are necessary to run this software:

- rdkit
- biopython
- openbabel

The following packages are also necessary in R (but should be autoinstalled the first time you run it):

- ggplot2
- SHAPforxgboost
- xgboost
- fastshap


#### Download the files ########################################################

Download all files from the DrugPred_RNA repository on GitHub and maintain the folder structure

#### Preparing the conda environment to run DrugPred_RNA ########################################################

It's possible to create a conda environment with the necessary packages from the .yml-file. To do this, run

conda env create --file DrugPred_RNA_env.yml

You should now be able to run

conda activate drugpred_RNA

#### Unzip the database containg drug molecules for docking  and prepare environment  ########################################################

1) create the environment variable DOCK_BASE_3_6 and point it to the folder containing DOCK 3.6 (e. g. setenv DOCK_BASE_3_6 /path_to_dock_folder/dockenv)

2) create the environment variable DOCK_BASE and point it to the same folder as $DOCK_BASE_3_6

3)  create the environment variable DrugPred and point it to the foler containing the DrugPred scripts

4) at the following folders to your path: $DOCK_BASE/bin/, $DOCK_BASE/private, $DOCK_BASE/etc, $DOCK_BASE/scripts, $DOCK_BASE/bin/Linux, $path.

(e. g. set path = ($DOCK_BASE/bin/ $DOCK_BASE/private $DOCK_BASE/etc $DOCK_BASE/scripts $DOCK_BASE/bin/Linux $path))


#### Downloading and preparing files for DrugPred_RNA ###########################################################

STEP 1.1:
cd into the main directoy of DrugPred_RNA. Here, you should have the folders "scripts" and "input_csvs".

STEP 1.2:
Create a .txt file, and name it "PDB_ids.txt". Add the four-letter
PDB ids you want DrugPred_RNA to predict in this document, separating entries with a comma. ie.
"2yie, 5c45" etc.

STEP 1.3:
Run:
python scripts/retrieve_pdb_structures.py -i PDB_ids.txt

Step 1.4:

cd into the directory named "PDB-FILES", which contains all the downloaded .cif-files with the matching PDB-id's.

Step 1.5:

Run:
python ../scripts/clean_up_pdb_files.py -f .cif

Step 1.6:

In the main directory there is now a file name "entries.csv", which contains the ids of all the entries which DrugPred_RNA will evaluate. The file consists of the PDB id and the three letter ligand code, separated by an underscore.
Note that if a binding site contains a metal ion, separate entries will be created for each individual metal, ie. "2yie_FMN", "2yie_FMN_MG". 

#### Docking, calculating and predicting the binding sites with DrugPred_RNA  ###################################

STEP 2.1:

cd into the directory PDB-FILES

STEP 2.2
Run:
python ../scripts/drugpred.py -s -1.3 --density 2

This may take a while as DrugPred will dock each binding site

#### Analyzing the output from DrugPred_RNA   ###################################################################

When the script is finished, cd back to the main directory. You should now have the following files:
 - "Predictions.csv", which contains the predicted outcomes of the assessed binding sites. Evaluated binding sites are evaluated as druggable (>= 0.5) or less druggable (<0.5)
 - "descriptor_values.csv", which contains the calculated descriptor valus for each assessed binding site.
 - In the "cleaned"-directory, subdirectories for each assessed binding site is made. Here, cd into the docking subdirectory to find:
	1. The binding site which was used for docking, with the suffix "_cofac.pdb"
	2. The original ligand of the binding site, with the suffix "_lig.pdb"
	3. The superligand file with the suffix "exp_atoms_modified.pdb".
	   This is the superligand created by DrugPred_RNA. 
   	   Atoms with "TRP"-labelled residues are solvent exposed, "ASP" are buried.
	4. The atoms which are in direct contact with the superligand has the suffix "_1.0_bs.pdb"
 - A directory named "ind_shap_plots". In here, SHAP plots for each assessed binding site is stored. These show the top 5 strongest influences descriptors have on the predicted outcome of each binding site.
