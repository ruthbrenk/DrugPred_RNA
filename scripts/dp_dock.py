#!/usr/bin/env python3
# should convert pdb file to format required for docking, prepares also co-factor files
#formerly known as setPDBatomsC.py
#I think we need Protein plus co-factor + metal, ligand file (AK had more, keep an eye on these changes)
#introduced check that ligand occurs only once in PDB file (RB)

import os,sys,string, os.path,math
import string,os, sys
#import MySQLdb
#from openeye.oechem import *
from argparse import ArgumentParser



#*******************************************************************************
# MAIN help
#**********


description = "Script to carry out the docking part of DrugPred_RNA"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

parser.add_argument("-prot",  action="store", dest="prot", help="PDB ID")
parser.add_argument("-id", "--id", type=str, action="store", dest="id", help="ID code of prot",required=True)
parser.add_argument("-v", "--dock_version",  action="store", dest="version", default='3.6', help="Dock version (3.5 or 3.6)",choices=['3.5', '3.6'])
parser.add_argument("-metal",  action="store", dest="metal",default=False, help="Metal in binding site", nargs='*')
parser.add_argument("-lig",  action="store", type = str, dest="lig",default=False, help="Ligand in binding site")
parser.add_argument("--density", action="store", dest="density", default='0', help="Superligand density, 0 = no check")
parser.add_argument("-s","--sl_cut_off_score", action="store", dest="sl_cut_off", default='-1.1', help="Score cut off for superligand, default = -1.1")
parser.add_argument("--descriptors", action="store", dest="des_type", default='new', help="old or new type of descriptors (modify SA for N/O)",choices=['new', 'old', 'test'])


args = parser.parse_args()
args_dict = vars(args)





#*******************************************************************************
# MAIN
#***************

#docking_scripts = os.environ['DockingScripts']
drugpred_scripts = "../../../scripts/"
#               "../../../../scripts/drug_moles.db"


#print args_dict
locals().update(args_dict) #generates local variables from key / value pairs

print("working on", id)
print(prot)
input_file_name = prot + '.pdb'
#check if file exits
if not os.path.isfile(input_file_name):
	#file does not exist
	print('file ', input_file_name, ' does not exist') 
	print('quitting dp_docking.py')
	sys.exit()
pdb_file = open(prot + '.pdb', 'r') #input file
out_file_protein= open(prot + '_C.pdb','w') #output file
out_file_protein_cofact= open(prot + '_cofac.pdb','w') #output file


ligand = lig
print(ligand, 'ligand')
out_file_ligand= open(prot + '_' + ligand + '_lig.pdb','w')

cofactor = " "
print(cofactor)

print('metal:', metal)
		


pdb_lines = pdb_file.readlines()

receptor = []
ligand_lines = []

for line in range(0,len(pdb_lines)):
	i = pdb_lines[line]
	find_H = i[34:].find('H') #search after RES NAME column
	if find_H == -1: #delete H atoms -> H not found
		i  = i[:20]+'   '+i[23:] #replace chain letter

	if i[:4]=="ATOM": #These lines contain the relevant information
		out_file_protein_cofact.write(i)
		outline="ATOM   "+i[7:11]+" "+ " CX" +"  "+"AAA"+""+i[20:54]+" \n"
		out_file_protein.write(outline)
		receptor.append(i)

	if i[:6]=="HETATM": #If structure contains cofactor, write it in pdb_output file
		if i[17:20].strip() == cofactor:
			i = str.replace(i, 'HETATM', 'ATOM  ') #just to be save
			out_file_protein_cofact.write(i)
			outline="ATOM   "+i[7:11]+" "+ " CX" +"  "+"AAA"+""+i[20:54]+" \n"
			out_file_protein.write(outline)
			receptor.append(i)
		elif i[17:20].strip() in metal:
			i = str.replace(i, 'HETATM', 'ATOM  ') #just to be save
			out_file_protein_cofact.write(i)
			outline="ATOM   "+i[7:11]+" "+ " CX" +"  "+"AAA"+""+i[20:54]+" \n"
			out_file_protein.write(outline)
			receptor.append(i)

		elif i[17:20].strip() == ligand:
			outline = str.replace(i, 'HETATM', 'ATOM  ') #HETATM must be ATOM in xtal-lig.pdb
			#print outline
			#out_file_ligand.write(outline)
			ligand_lines.append(outline)

print(len(receptor), "lines in receptor file")
print(len(ligand_lines), "lines in reduced sphere set")

if len(receptor)==0 or len(ligand_lines) ==0:
	print("len(receptor)==0 or len(ligand_lines) ==0")

	print('quitting dp_dock.py')
	cursor.close()
	conn.close()
	sys.exit()


#check distance
for atom in ligand_lines:
	x_lig = float(atom[30:37])
	y_lig = float(atom[38:45])
	z_lig = float(atom[46:53])
	#print x_lig, y_lig, z_lig, 'coords'
	#print atom
	for prot_atom in receptor:
		x_prot = float(prot_atom[30:37])
		y_prot = float(prot_atom[38:45])
		z_prot = float(prot_atom[46:53])
		dist = math.sqrt (math.pow((x_lig-x_prot),2) + math.pow((y_lig-y_prot),2) + math.pow((z_lig-z_prot),2))
		#print dist, 'dist'
		if dist < 5.0:
			out_file_ligand.write(atom)
			break #no need to check for more distances, atom is close enough

out_file_ligand.close()
out_file_protein.close()
out_file_protein_cofact.close()


command = 'cp ' + prot + '_' + ligand + '_lig.pdb' + ' xtal-lig.pdb'
print(command)
os.system(command)

#command = 'touch xtal-lig.mol2' 
#os.system(command)

		
#create rec.amb
command = "python3 " + drugpred_scripts + 'pdb2amb.py ' + prot + '_C.pdb rec.amb'
print(command)
os.system(command)

#create link from rec.amb to rec.pdb
command = 'ln -s rec.amb rec.pdb'
print(command)
os.system(command)

#make file

print(' make file:')
command = 'make auto -f ' + drugpred_scripts + '/Makefile_DrugPred_' + version
print(command)
os.system(command)

#check if make file went through

if not os.path.isdir('testing'): #make file failed
	print('Makefile failed')
	print('STOP!!!!!!!')
	sys.exit()

#check how many sphes, if at least 40 => use INDOCK file with faster settings

#sph_file = open('sph/match.sph', 'r')
#sph_lines = sph_file.readlines()
#num_sph = len(sph_lines)

#if num_sph < 41: #first line in match.sph is header
	#indock_file = 'INDOCK'
#else:
	#indock_file = 'INDOCK_fast'
#sph_file.close()

indock_file = 'INDOCK'
#INDOCK

command = 'cp ' + drugpred_scripts + '/' + indock_file + ' testing'
print(command)
os.system(command)


#do docking
os.chdir('testing')
path = "../" + drugpred_scripts + '/drug_moles.db'
if os.path.isfile(path):
	print("Drug Moles exist")
else:
	print("Wrong path")

command = 'ln -s ' + "../" + drugpred_scripts + '/drug_moles.db .'

os.system(command)

#Mark that docking has started

print('docking started')

command = os.environ['DOCK_BASE' + '_' + version.replace('.','_')] + 'bin/Linux/dock.csh'

print(command)

os.system(command)

outdock_file = open('OUTDOCK', 'r')
outdock_content = outdock_file.read()
if outdock_content.find('too many points in force field grid') != -1:
	print("Docking failed")


os.chdir('../') #back to directory docking

print('finished docking')

