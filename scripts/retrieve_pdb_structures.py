"""
This piece of code uses BioPython to download .pdb-files from a list of pdb-codes ("RNA_w_lig_PDB_codes") and sends them to a directory labeled with today's date.
It also identifies the ligand (3-letter code) and metals that are within 5Å of the ligand, and uploads these descriptors
to a MySQL-table.

Common buffer components (such as DMSO, glycol, etc.) are filtered out by passing the heteroresidues through a common buffer
The distinction between ligands and modified residues are tricky, as both are classified as heteroligands by pdb-files and BioPython.
Entries like 6fz0 can be tricky, as the ligand (SAM) here also appears as a modified residue, and would normally be kicked out.
This script has a two-fold approach to dealing with this problem. 

If A) The .pdb-file contains a MODRES-line:
	the script creates a list of these known modified residues, and any heteroresidue which is not in this list and also not a common buffer component is defined as a ligand

or B) The .pdb does not contain a MODRES-line:
	the script filters out all the heteroresidues which are both in the common buffer component list and the common modified residue list

by Illimar Rekand,
Brenk group
October 2019
illimar.rekand@uib.no
illimar.rekand@gmail.com

"""

import Bio.PDB
import os, sys, urllib, string, tarfile
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqIO import PdbIO
from datetime import date
from argparse import ArgumentParser
import csv
import re      #for filtering out words

#print("BIOPYTHON VERSION", Bio.__version__)


dir_path = os.path.dirname(os.path.realpath(__file__))
local_dir = "./PDB-FILES"
cleaned = "./PDB-FILES/cleaned"
if os.path.exists(local_dir):
	print(local_dir + " already exists!")
else:
	mk_local_dir = os.mkdir(local_dir)
	mk_cleaned_dir = os.mkdir(cleaned)
#import my_mysql3 as mysql

description = "Retrieve PDB files and create a list of input lines"
parser = ArgumentParser(description=description)
parser.add_argument("-i", "--input_file", type=str, action="store", dest="input_list", help="Name of csv file with PDB files",required=True)
#parser.add_argument("-t", "--data_table", type=str, action="store", dest="tb",  help="table which contains input DrugPred calculations",required=True)
#parser.add_argument("-u", "--user", type=str, action="store", dest="us", help="MySQL usernmae",required=True)
#parser.add_argument("-p", "--password", type=str, action="store", dest="pw", help="MySQL password",required=True)
#parser.add_argument("-f",  type = str, dest="format", help="run on cluster")#

args = parser.parse_args()
args_dict = vars(args)
print (args_dict)
locals().update(args_dict) #generates local variables from key / value pairs

entry_id_list = []
idx = 0

#-----------------------------------
def write(entry_id, metal, residues):
	global entry_id_list
	global idx
	entry_id_list.append(entry_id)
	idx += 1
	string = str(idx) + " \t" + entry_id + "\t" + str(structure_ID).lower() + "\t" + str(residues) + "\t" + str(metal) + "\n"
	print(string)
	out.write(string)

def get_codes(csv_name):
	#get  3 letter salt_codes for metals, buffer components, etc.
	code_list = []
	with open("./input_csvs/" + csv_name) as csvfile:
		list_from_csv =  [row[0] for row in csv.reader(csvfile, delimiter= ",")]
	return list_from_csv

def update_entry_row(entry_id):
	command = "Select * from "+ tb +" where id = '" + entry_id  + "'"  #Check if entry in MySQL already exists
	print(command)
	cursor.execute(command)
	label_exists = cursor.fetchall()
	if len(label_exists) > 0: #Delete entry if this already exists, so we can update it
		print(label_exists[0][0], label_exists[0][1])
		command = "delete from " + tb + " WHERE id = '" + entry_id + "'"
		print(command)
		cursor.execute(command)

def find_metals(model):
	for chain in model:
		for residue in chain:
			residue_id = residue.get_id()
			hetfield = residue_id[0]
			if hetfield[0]=="H": #this means the residue is a heteroresidue
				res_name = residue_id[0][2:].strip()
				for atom in residue:
					if (atom.get_name() in salt_codes) and atom.get_name() not in ["S", "O1", "O2", "O3", "O4", "CL" ]: #this atom is probably a metal-atom, later we will check if it's in the binding site
						metal_atom_list.append(atom)

def inspect_chains_and_upload(model):
	global idx
	for chain in model:
		print(metal_atom_list)
		residue_list=[]
		chain_string = str(chain)
		for residue in chain:
			residue_id = residue.get_id()
			hetfield = residue_id[0]
			if hetfield[0]=="H":
				res_name = residue_id[0][2:].strip()
				#if res_name.strip() not in buffer_codes and res_name.strip() not in salt_codes and res_name.strip() not in modres_list and res_name.strip() not in residue_list:
				if res_name.strip() not in (buffer_codes + salt_codes + modres_list + residue_list): #this residue is neither a buffer-component, a metal, a modified residue in the polymer chain or already accounted for
					print(res_name)
					residue_list.append(res_name)
					if len(metal_atom_list) > 0: 
						atom_list = []
						for atom in residue:
					  		atom_list.append(atom)
						if len(atom_list) > 0:
							for metal_atom in metal_atom_list:
								for atom in atom_list:
									distance = atom - metal_atom
									if distance < 5: #is metal atom 5Å away from ligand?
										bs_metals.append(metal_atom.get_name())
										break #metal is close to this residue, no need to check other atoms in this resi

					if len(bs_metals) == 0:
						print("No metals found")


		print(structure_ID)
		print(bs_metals)
		print(residue_list)


		command_list = []
		

		first_entry = False

		for residues in set(residue_list):
			if len(bs_metals) > 0: #First makes a "normal" entry, then adds subsequent metal entries. Elegant? No. Works? Yes.
				entry_id = str(structure_ID).lower() + "_" + str(residues) #the entry id-for the non-metal entry
				#update_entry_row(entry_id)
				#command = "INSERT INTO " + tb + " (id, prot, lig, resolution, comment) values ('" + entry_id + "', '" + str(structure_ID).lower()  + "', '" + str(residues) + "', '" + str(resolution) + "', '" + name + "' )" #make one entry in MySQL without metal
				if entry_id not in entry_id_list:
					write(entry_id, " ", residues)
				for metal in set(bs_metals): #set() prevents multiple entries of same metal, this would raise MySQL-error 1062
					entry_id = str(structure_ID).lower()  + "_" + str(residues) + "_" + str(metal) #the entry id for metal-containing structure
					print(entry_id)
					if entry_id not in entry_id_list:
						write(entry_id, metal, residues)

					if structure_ID not in Unique_structures:
						Unique_structures.append(structure_ID)
			elif len(bs_metals) == 0:
				entry_id = str(structure_ID).lower()  + "_" + str(residues)
				#update_entry_row(entry_id)
				#command = "INSERT INTO " + tb + " (id, prot, lig, resolution, comment) values ('" + entry_id + "', '" + str(structure_ID).lower()  + "', '" + str(residues) + "', '" + str(resolution) + "', '" + name + "' )" #make one entry in MySQL without metal
				if entry_id not in entry_id_list:
					write(entry_id, " ", residues)
				if structure_ID not in Unique_structures:
					Unique_structures.append(structure_ID)
#-----------------------------------


f1= open(input_list,"r") #Opens file with PDBs

for line in f1:
	print(line, type(line))

PDB_list = line.split(", ") #Makes a list of the PDB-files input .txt-file
print(PDB_list)
last_entry = PDB_list[-1]
if len(last_entry) > 4: #fixes the /n in the last entry in the PDB-list. Annoying.
	last_entry_fixed = last_entry[0:4]
	PDB_list[-1] = last_entry_fixed

#print mod_bases_read, "------------------------------------------------------"

#print mod_bases_list, "******************************************************"

#DB connection
#conn=mysql.connect2server('', 'webuser','drug_pred')
#cursor = conn.cursor ()  	


#get buffer component code
buffer_codes = get_codes('crystal_buffer_components.csv')   # Makes list of three-letter buffer codes, in the drug_pred db
salt_codes = get_codes('metal_codes.csv')                   # makes list of metal codes, in the drug_pred db
print(buffer_codes)
print(salt_codes)
SO4 = ["S", "O1", "O2", "O3", "O4" ] #apparently a lot of entries (2oe5, )
salt_codes = salt_codes + SO4
mmcif_failed_list = []
Unique_structures = []


#print(commonhetres_list)
#----------------------------------

out = open("entries.csv", "w+")
with open("./input_csvs/mod_res.csv") as csvfile:
		commonhetres_list =  [row[0] for row in csv.reader(csvfile, delimiter= ",")]

for structure_ID in PDB_list:
	#db = str(sys.argv[2]) #ie. drugpred_rna_ire
	#tb = str(sys.argv[3]) #the MySQL-table that accepts the pdb-codes, lig-codes and metals. ie. RNA_structures_3, test_input_files
	#us = str(sys.argv[4]) # your username
	#pw = str(sys.argv[5]) # your password
	#print(db, us, pw)


	#conn=mysql.connect2server(pw, us, db)  			    
	#cursor = conn.cursor()

	#commonhetres_list = get_codes("codes", "drugpred_rna_ire.modified_bases_rna_dna") #Fetches the common modified residues list, in the drug_pred_rna_ire db

	structure_id = structure_ID.lower()    # Making pdb-code from upper to lower case
	#print(structure_id)m
	 
	####Retrieve pdb####
	pdbl = PDBList()
	filename_list = []
	#filename = pdbl.retrieve_pdb_file(structure_id, pdir=local_dir, file_format="pdb") #retrieves the PDB-file, saves it in a directory with appropiate date. File is saved with pdb
	filename = pdbl.retrieve_pdb_file(structure_id, pdir=local_dir, file_format="mmCif") #retrieves the equivalent mmCif-file, we need this to look at polymerized vs nonpolymerized heteroresidues
	filename_list.append(filename)
	if os.path.exists(filename):
		print(filename + " does exist")
	else:
		print(filename + " does not exist")
		filename = pdbl.retrieve_pdb_file(structure_id, pdir=local_dir, file_format="bundle")
		print(filename, "filename")
		count = 0
		with tarfile.open(str(filename)) as tf: #open tarfile
			count = sum(1 for member in tf if member.isreg()) #count every file inside the tar
		tf = tarfile.open(str(filename), "r:")
		tf.extractall(path=local_dir)
		tf.close()
		count -= 1
		print("count", count)
		filename_list = [local_dir + "/" + structure_id + "-pdb-bundle" + str(x + 1) + ".pdb" for x in range(count)] #create a list of the directories for each pdb inside the tar-bundle
		print("filename_list", filename_list)
	print("filename", filename)
	filename = pdbl.retrieve_pdb_file(structure_id, pdir=local_dir, file_format="mmCif") #retrieves the equivalent mmCif-file, we need this to look at polymerized vs nonpolymerized heteroresidues

	print(filename)
		#---------------parse the header-------------------------------------------------
	parser = PDBParser(PERMISSIVE=1)
	parser_cif = MMCIFParser() 
	structure = None
	multiple_PDBs = 0
	try:
		structure_cif = parser_cif.get_structure(structure_id, filename)
	except:
		print("MMCIF parsing failed")
		mmcif_failed_list.append(structure_ID)
		continue


	mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
	sc = mmcif_dict['_entity_poly.nstd_monomer']

	mon_ids = mmcif_dict['_entity_poly_seq.mon_id']
	hetero = mmcif_dict['_entity_poly_seq.hetero']

	modres_list = [] #This is a list of the known modified residues that reside within a polymer chain inside this speicific structure

	for i in range(0,len(mon_ids)):
		if (mon_ids[i] in commonhetres_list) and (mon_ids[i] not in modres_list):
			modres_list.append(mon_ids[i])
	print("modres_list", modres_list)


	#---------------get info from .pdb file-------------------------------------------------
	for filename in filename_list:
		structure = parser_cif.get_structure(structure_id, filename)
		print("filename", filename, structure)
		#resolution = structure.header["resolution"]
		#deposition_date = structure.header["deposition_date"]
		#compound = structure.header["compound"]
		#name = structure.header["name"]
		model = structure[0] #only takes the first model in the entry, NMR have several entries. The cleanup script also only cleans the first model
		bs_metals = []
		metal_atom_list=[] #list of all metals in the pdb-file
		find_metals(model)
		inspect_chains_and_upload(model)


print(Unique_structures)
print(entry_id_list)
print("Number of unique structures:", len(Unique_structures))
print("Structures where mmCif-file failed", mmcif_failed_list)

out.close()
#-----------------------------



