#!/usr/bin/env python3

#master script to run DrugPred
#cmd line version


#setup directory for each drugpred job
#prepare for docking
#dock
#receptors: must be called 'prot'.pdb


import sys as os, os.path
from argparse import ArgumentParser
import csv

#*******************************************************************************
# MAIN help
#**********
#


description = "Script to run DrugPred 2.0"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

#parser.add_argument("-d", "--database", type=str, action="store", dest="db", help="database which contains input and output table for DrugPred",required=True)
#parser.add_argument("-t", "--data_table", type=str, action="store", dest="tb",  help="table which contains input DrugPred calculations",required=True)
#parser.add_argument("-r", "--result_table", type=str, action="store", dest="dt",  help="table which contains results of DrugPred calculations",required=True)
#parser.add_argument("-u", "--user", type=str, action="store", dest="us", help="MySQL usernmae",required=True)
#parser.add_argument("-p", "--password", type=str, action="store", dest="pw", help="MySQL password",required=True)
parser.add_argument("-c", "--cluster",  action="store_true", dest="cluster", help="run on cluster")
parser.add_argument("-m", "--method",  type=str, action="store", dest="method", default='all', help="What shall I do?",choices=['all', 'docking', 'superligand', 'descriptors','overlapp', 'predict'])
parser.add_argument("-v", "--dock_version",  action="store", dest="version", default='3.6', help="Dock version (3.5 or 3.6)",choices=['3.5', '3.6'])
parser.add_argument("--delete",  action="store_true", dest="delete",default=False, help="Delete directory if output already exists?")
parser.add_argument("--density", action="store", dest="density", default='0', help="Superligand density, 0 = no check")
parser.add_argument("-s","--sl_cut_off_score", action="store", dest="sl_cut_off", default='-1.1', help="Score cut off for superligand, default = -1.1")
parser.add_argument("--debug",  action="store_true", dest="debug",default=False, help="Print debugging output")

args = parser.parse_args()
args_dict = vars(args)
locals().update(args_dict) #generates local variables from key / value pairs


def predict():
	if not os.path.isdir("../ind_shap_plots"):
		os.mkdir("../ind_shap_plots")
	command = "Rscript ../scripts/predict.r"
	print(command)
	os.system(command) 


#*******************************************************************************
# MAIN
#***************

drugpred_path = "../scripts/"
print(drugpred_path)
id_list = []




#print (args_dict)

print (version)

#print delete

#check that db values are okay
#tables
#conn=mysql.connect2server(pw, us, db)  			    
#cursor = conn.cursor ()

#check environment variables


#---------------------------------
#---------------------------------

#get to do from table


with open("../entries.csv") as csvfile:
	entries = csv.reader(csvfile, delimiter= "\t")
	entries = [tuple(row) for row in entries]

print ('delete existing files:', delete, 'run method: ', method)
counter = 0

count_thousands = 0
new_number = 1
old_number = 1
size = 1000
count_dict = {}

with open("../descriptor_values.csv", "w+") as out:
	close = True
	writer = csv.writer(out, delimiter = ",")
	#out.write("id \t csa \t sl_csa \t hsa \t exp_sl_sa \t fsasa \t asa \t asa_r \t chsa \t chsa_r \t psa_r \t aliphat_aromat_t \t fr_hpb_atoms \t c_ali_sa_r \t ali_sa_r \t  binding_site_seq \t  SpherocityIndex \t  Asphericity \t  Eccentricity \t  InertialShapeFactor \t  RadiusOfGyration \t  NPR2 \t  NPR1 \t  PMI3 \t  PMI2 \t  PM1 \t  no_sl_atoms \t  no_bs_atoms \t  fr_buried_sl_atoms \t  sl_bs_r \t  fr_surf_sl_atoms \t  vol \t  csa_vol_r \t  sa_vol_r \n")
	writer.writerow( ["id" ,   "csa" ,   "sl_csa" ,   "hsa" ,   "exp_sl_sa" ,   "fsasa" ,  "chsa" ,   "chsa_r" ,   "psa_r" ,   "aliphat_aromat_t" ,   "fr_hpb_atoms" ,   "c_ali_sa_r" ,   "ali_sa_r" ,   "binding_site_seq" ,   "SpherocityIndex" ,   "Asphericity" ,   "Eccentricity" ,   "InertialShapeFactor" ,   "RadiusOfGyration" ,   "NPR2" ,   "NPR1" ,   "PMI3" ,   "PMI2" ,   "PM1" ,   "no_sl_atoms" ,   "no_bs_atoms" ,   "fr_buried_sl_atoms" ,   "sl_bs_r" ,   "fr_surf_sl_atoms" ,   "vol" ,   "csa_vol_r" ,   "sa_vol_r" ] )

if method == 'predict':
	if os.path.isfile('../descriptor_values.csv'):
		predict()
		print("Predictions finished")
		exit()
	else:
		print("Descriptor value file does not exist")

for idx, id, prot, lig, metal in entries:
	if id not in id_list:
		id_list.append(id)
	else:
		continue
	print(id, 'working with this one')
	#check if pdb file exists
	if not os.path.isfile(prot + '.pdb'):
		print (prot, ' prot does not exist')
		continue
	if os.path.isdir(id): 
		if delete and method == 'all': #delete old directory
			os.system('rm -rf ' + id)
		else:
			print (id, ' does already exist')
			if not delete and method == 'all':
				print ('skip')
				continue
	elif method != 'all':	#can not caluculate des if directory is not there
		if method != 'docking':
			print (id, ' id does not exist')
			continue
	
	print(cluster, 'cluster')
	if not cluster:
		print('running locally')
		print(metal)
		command = 'python3 ' + drugpred_path + 'drugpred_call_functions.py -id ' + str(id) +  ' -prot ' + str(prot) + ' --dock_version ' + version + ' -m ' + str(method) + ' --density ' + density + ' -s ' + sl_cut_off + ' -metal ' + metal + ' -lig ' + lig 
		print (command)
		os.system(command)
	else: #cluster

		new_number = count_thousands / size + 1 #=1 necessary, otherwise first file will always be called 0
		if new_number != old_number :
			old_number = new_number
			count_dict[new_number-1] = counter
			counter = 0

		counter = counter + 1
		count_thousands = count_thousands + 1

		#print 'submit job on cluster'

		first_bit = "#$ -S /bin/tcsh\n#$ -cwd\nworkdir=/$SCRATCH/$SLURM_ARRAY_TASK_ID\nmkdir -p $workdir\ncd $SLURM_SUBMIT_DIR\n"

		file_name = str(counter) + '_' + str(old_number) + '_start.bin'
		start_file = open(file_name, 'w')
		start_file.write(first_bit)

		start_file.write('echo $workdir\n')

		start_file.write('cp ' + prot + '.pdb $workdir\n')

		# print os.path.isdir(id)
		if os.path.isdir(id): #dir got already delete above if necessary
	        	start_file.write('cp -R ' + id + ' $workdir\n')

		start_file.write('cd $workdir\n')
		start_file.write('$DrugPred/drugpred_call_functions.py -db ' + db + ' -tb ' + tb + ' -dt ' + dt  + ' -user ' + us + ' -password '+ pw + ' -id ' + id +  ' -prot ' + prot +  ' -dock_version ' + version + ' -m ' + method +  ' --density ' + density + ' -s ' + sl_cut_off + ' -des_type ' + des_type + ' -debug ' + str(debug) + '\n')
		start_file.write('cp -R ' + id  + ' $SLURM_SUBMIT_DIR\n')
		
		if not debug:

	                start_file.write('rm -rf  $SLURM_SUBMIT_DIR/' + id + '/docking/sph\n')
        	        start_file.write('rm -rf  $SLURM_SUBMIT_DIR/' + id + '/docking/grids')

		start_file.close()

		os.system('chmod 744 '  + file_name)

if method in ('all'):
	predict()



if cluster: #submit job
	count_dict[new_number] = counter
	print (count_dict)
	for i in range(1,new_number+1): 

		file_name = str(i) + '_' + 'array_start.bin'
		start_file = open(file_name, 'w')
		start_file.write("#!/bin/bash\n")
		start_file.write('#SBATCH --time=0-03:00:00\n')
		start_file.write('#SBATCH --partition normal\n')
		start_file.write('${SLURM_ARRAY_TASK_ID}_' + str(i) + '_start.bin\n')
		start_file.close()
		os.system('chmod 744 '  + file_name)

		print (i)
		command = 'sbatch --array=1-' + str(count_dict[i]) + ' ' + file_name  
		print (command)
		os.system(command)


if close:
	out.close()
print ('done!!!!!!!!!!!!!')

#to do:
#check with bigger super ligand if I get about the same descriptors
#optimize docking setup
#make ready for cluster
#redock test set
#dock RNA molecules
#analyse
#write paper








