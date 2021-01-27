#!/usr/bin/env python3

#run docking calculations, superligand and calculate descriptors

import sys, os, os.path, csv
from argparse import ArgumentParser

#drugpred_call_functions.py -db [data base] -tb [data table] -dt [result table] -user [username db] -password [password] -id [id] -prot [prot] -dock_version [version]"

drugpred_path = "../../../scripts/"

description = "Script to run DrugPred_RNA"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

parser.add_argument("-id",  action="store", dest="id", help="Entry id")
parser.add_argument("-prot",  action="store", dest="prot", help="PDB ID")
parser.add_argument("-metal",  action="store", type = str, dest="metal", nargs='*', default=False, help="Metal in binding site")
parser.add_argument("-lig",  action="store", type = str, dest="lig",default=False, help="Ligand in binding site")
parser.add_argument("-c", "--cluster",  action="store_true", dest="cluster", help="run on cluster")
parser.add_argument("-m", "--method",  type=str, action="store", dest="method", default='all', help="What shall I do?",choices=['all', 'docking', 'superligand', 'descriptors','overlapp', 'predict'])
parser.add_argument("-v", "--dock_version",  action="store", dest="version", default='3.6', help="Dock version (3.5 or 3.6)",choices=['3.5', '3.6'])
parser.add_argument("--delete",  action="store_true", dest="delete",default=False, help="Delete directory if output already exists?")
parser.add_argument("--density", action="store", dest="density", default='0', help="Superligand density, 0 = no check")
parser.add_argument("-s","--sl_cut_off_score", action="store", dest="sl_cut_off", default='-1.1', help="Score cut off for superligand, default = -1.1")
parser.add_argument("--debug",  action="store_true", dest="debug",default=False, help="Print debugging output")


args = parser.parse_args()
args_dict = vars(args)

print(args_dict)
locals().update(args_dict) #generates local variables from key / value pairs

if metal:
	metal = metal[0]
else:
	metal = 'XX'

#---------------------------------
def docking(id,prot):
	if not os.path.isdir(id):
		os.mkdir(id)
		os.chdir(id)
		os.mkdir('docking')
	else:
		os.chdir(id)

	#write out version
	tmp_file =open('DOCK_' + version + ".txt",'w')
	tmp_file.close()

	os.chdir('docking')



	command = 'cp ../../' + prot + '.pdb' + ' .'
	print(command)
	os.system(command)
	#this script does the docking stuff
	command = 'python3 ' + drugpred_path + 'dp_dock.py -id ' + str(id) +  ' -prot ' + str(prot) + ' --dock_version ' + version + ' -m ' + str(method) + ' --density ' + density + ' -s ' + sl_cut_off + ' -metal ' + str(metal) + ' -lig ' + lig
	print(command)
	os.system(command)
	os.chdir('../..')

#---------------------------------
def superligand(id, density,sl_cut_off):
	os.chdir(id)
	os.chdir('docking')
	os.chdir('testing')
	if os.path.isfile('res.eel1.gz'): #docking produced an output
		os.system('gunzip res.eel1.gz')
	else: #nothing was there, generate empty file
		os.system('touch res.eel1')
	#calculate superligand
	command = "python3 ../" +drugpred_path + 'superligand.py ' + density + ' ' + sl_cut_off
	print(command)
	os.system(command)
	os.chdir('../../..')		


#---------------------------------
def descriptors(id):

	os.chdir(id)
	os.chdir('docking')

	command = "python3 " + drugpred_path + 'calculate_descriptors.py -id ' + id + ' -prot ' + prot 

	
	print(command)
	print('command', command)
	os.system(command)
	os.chdir('../..')
#---------------------------------
def overlapp(id):
	os.chdir(id)
	os.chdir('docking')
	os.chdir('testing')
	command = drugpred_path + 'overlapp_superligand_xtal_lig.py -db ' + db + ' -tb ' + tb + ' -user ' + us + ' -password '+ pw + ' -id ' + id + ' -dt ' + dt
	print(command)
	os.system(command)
	os.chdir('../../..')		
#---------------------------------

if os.path.exists(id) and method =='all':
	print('dir', id, ' does already exist, continue with next receptor')
	sys.exit()

if method == 'all' or method == 'docking':

	#prepare for docking and dock
	docking(id,prot)
if method == 'all' or method == 'superligand':
	#Generate superligand
	superligand(id,density,sl_cut_off)

if method == 'all' or method == 'descriptors':

	#calculcate descriptors
	descriptors(id)

if method == 'overlapp':

	#calculcate overlapp
	overlapp(id)


print('finished calculations!!!!!!!!!!')

