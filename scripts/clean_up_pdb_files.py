#!/usr/bin/env python3

"""
This script carves out the binding site around a pre-determined ligand, which is fetched from a MySQL-database. It also removes water molecules, metals and so-called buffer components.
It inputs a downloaded pdb-file from a table (e.g. test_input_files or RNA_structures), and outputs a pdb-file. The ligands should be outputted as HETATM lines in the final pdb-file.
This script also saves only the first conformer of the pdb-file, which is crucial for RDKit calculcations later on.

KNOWN BUGS and issues:
* a lot of pdb-files are inconsistently formatted. There has been a lot of trouble with whitespaces in strings.

Illimar Rekand,
November 2019

illimar.rekand@uib.no
illimar.rekand@gmail.com


"""
import os, sys, string


from Bio.PDB import *
from Bio import PDB
from argparse import ArgumentParser
import csv

import numpy as np
#import xpdb

global keep_res_list


clean_cofac = False #outputs the cleaned file in the "docking" subdirectory
send_to_clean_dir = True #Outputs the cleaned file to a "cleaned" subdirectory. This 
testing = False
get_first_model = True #gets the first model of a pdb-file. Useful for NMR structures
read_in_pdb_files = True

#------------------------
class save_res(Select):
	def accept_residue(self, residue):
	#print '.................'
	#print residue.id
		#print(keep_res_list)
		if residue in keep_res_list:
			#print(residue)
			return 1
		else:
			return 0

#------------------------
description = "Clean up .pdb-files or .cif files to prepare them for docking"
parser = ArgumentParser(description=description)
parser.add_argument("-f",  type = str, dest="format", help="file format")

args = parser.parse_args()
args_dict = vars(args)
print (args_dict)
locals().update(args_dict) #generates local variables from key / value pairs


input_dir = "./"
input_file_format = format
if input_file_format in ["cif", "mmcif", ".cif", ".mmcif"]:
	read_in_pdb_files = False
elif input_file_format in ["pdb", ".pdb"]:
	read_in_pdb_files = True
else:
	print("Format undefined")
	sys.exit()


files = os.listdir(input_dir)


#conn=mysql.connect2server(password, username,'drugpred_rna_ire')
#cursor = conn.cursor()


#get list with modified bases

input_csv_path = "../input_csvs/"

with open(input_csv_path + "mod_res.csv") as csvfile:
	mod_res_list =  [row[0] for row in csv.reader(csvfile, delimiter= ",")]

#print (mod_res_list)

#get metals

with open(input_csv_path + "metal_codes.csv") as csvfile:
	metal_list =  [row[0] for row in csv.reader(csvfile, delimiter= ",")]
	#print(metal_list)


with open(input_csv_path + "crystal_buffer_components.csv") as csvfile:
	crystal_buffer_components_list =  [row[0] for row in csv.reader(csvfile, delimiter= ",")]
	#print(crystal_buffer_components_list)

print(crystal_buffer_components_list)
files_cleaned = 0
failed_save = []

crystal_pdb_list = []

for file in files:
	#if file[-3:] == 'cif':
	if file[-3:] == 'pdb' or file[-3:] == 'cif':
		print(file[0:4])
		code =file[:-4]
		

		with open("../entries.csv") as csvfile:

			lig_list =  [row[3] for row in csv.reader(csvfile, delimiter= "\t")]

		bool_list = [False for x in range(len(lig_list))]

		lig_dict = dict(zip(lig_list, bool_list)) # this way we ensure that the lig is properly stripped for if-statement later on

		if len(lig_dict) == 0: #nothing found, next pdb code
			continue

		keep_res_list = []
		used_chain_ids = {}

		if read_in_pdb_files:
			parser = PDBParser()
		else:
			parser = PDB.MMCIFParser() 
		print(input_dir + file)
		structure = parser.get_structure('xxx', input_dir + file)
		print("structure", structure)

		if get_first_model: #only get the first model of a pdb-file. Crucial for NMR files, otherwise all models will be combined into one structure later on
			#print (len(structure), "length structure")
			structure = structure[0]
			#print("Get first model == True")

		rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0])) # Reference frame for translating the coordinates closer to origo. If xyz-coordinates too large in pdb file, docking fails


		for atom in structure.get_atoms():
			atom_C1 = atom.coord.copy()
			break


		atom_list = Selection.unfold_entities(structure, 'A') # A for atoms
		ns = NeighborSearch(atom_list)
		# Iterate over all residues in a model

		for residue in structure.get_residues():
			#print("residue.get_resname", code, residue.get_resname(), residue.get_resname().strip())
			res_name = residue.get_resname().strip()
			if res_name in lig_dict: #only keep first instance of a ligand
				#print(lig_dict)
				if not lig_dict[residue.get_resname().strip()] : #first time

					lig_dict[residue.get_resname()] = True #ligand already found, don't keep again
					chain_id = residue.get_parent()
					print("PARENT CHAIN", chain_id)
					keep_res_list.append(residue)


					#print (keep_res_list, lig_dict)
					#print ("residue name", residue.get_resname())
					#print ('here')
			
					for atom in residue:
						#print(atom, residue)
						atom_line = atom.get_full_id()
						#print("atom_line", atom_line)
						center = atom.get_coord()
						neighbors = ns.search(center, 15.0) # 15.0 for distance in angstrom (electrostatic is not considered during docking, do not keep too much)
						residue_list = Selection.unfold_entities(neighbors, 'R') # R for residues
						for res in residue_list:
							if (res not in keep_res_list) and ((res.get_resname() not in lig_dict.keys()) or ( res_name not in lig_dict.keys()) or (res_name not in crystal_buffer_components_list)) and (res.get_resname() != 'HOH') : #exclude also water molecules
								#print("res", res, res_name, res.get_resname())
								if res.get_resname() not in crystal_buffer_components_list:
									keep_res_list.append(res)
								else:
									if code not in crystal_pdb_list:
										crystal_pdb_list.append(code)

			for atom in residue:
				atom.transform(rotation_matrix, -atom_C1) #every atom in structure moved closer to origo

		io = PDBIO()
		io.set_structure(structure)

		
		for keep_res in keep_res_list:
			print("kept res", keep_res.get_resname())
		try:
			io.save('out.pdb', save_res()) #creates an output file, which we redirect and edit down below. This file is overwritten for every file that is cleaned
		except:
			failed_save.append(code)
			continue

		#io.save('out.pdb')

		#change hetatm -> for linked res

		else:
			output_dir = "."
		pdb_file = open('out.pdb', 'r')
		pdb_out = open(output_dir + "/" + code +'.pdb','w')
		coord_list = []
		line_number = 1
		for line in pdb_file.readlines():
			skip = False
			#print("line[17:20]", line[17:20].strip(), "line[0:6]", line[0:6])
			if line[17:20].strip() in lig_dict.keys() and line[0:6] == "HETATM":
				if line[16] != " ": #removes the A, B or other conformer tag if the ligand has multiple conformers
					line = line[0:16] + " " + line[17:]
				line_number_string = str(line_number)
				line = line[0:6] + line_number_string.rjust(5) + line[11:] #rewrites the line number (fields 7-11) to make sure that there are no gaps. OpenBabel uses the line to index the atoms, while RDKit probably counts the number of lines in the pdb-file
				pdb_out.write(line)
				line_number += 1
				skip = True #doesn't write the ligand twice
			if line[17:20].strip() in metal_list and line[0:6] == "HETATM":
				print("Found metal")
				if line[16] != " ": #removes the A, B or other conformer tag if the ligand has multiple conformers
					line = line[0:16] + " " + line[17:]
				line_number_string = str(line_number)
				line = line[0:6] + line_number_string.rjust(5) + line[11:] #rewrites the line number (fields 7-11) to make sure that there are no gaps. OpenBabel uses the line to index the atoms, while RDKit probably counts the number of lines in the pdb-file
				pdb_out.write(line)
				line_number += 1
				skip = True #doesn't write the ligand twice
			#print(line[12:14])
			if not line[12:14].strip():
				skip = True
			if (line[17:20].strip() in mod_res_list) and (line[17:20].strip() not in lig_dict.keys()): #change modfied res (but only if this is not the current ligand)
				line = 'ATOM  ' + line[6:]
			if len(line) > 15: #Prevents same atom in two conformer outside the polymer chain to be written twice. Example 3e5e. Docking will work, but descriptors script will fail
				if line[16] != " ":
					if line[30:53] in coord_list:
						skip = True
					else:
						coord_list.append(line[30:53]) 

			if skip == False:
				line_number_string = str(line_number)
				line = line[0:6] + line_number_string.rjust(5) + line[11:] #same as above
				pdb_out.write('ATOM  ' + line[6:])
				line_number += 1
		pdb_file.close()
		pdb_out.close()
		files_cleaned += 1



print("Number of files cleaned:", files_cleaned)
print("Number of failed saves", failed_save)

os.system('rm -rf out.pdb') #removes the output-file, which is overwritten every time a new file is cleaned

