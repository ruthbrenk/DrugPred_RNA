#!/usr/bin/env python3

#read OUTDOCK file
# keep codes that fullfil score cut-off
#generate pdb file for these molecules
#thin out superligand

import string, os,math,sys


atom_cut_off = 1.2 #C-C distance
min_hits = int(sys.argv[1]) #how many atoms have to be found close to one atom

#score_cut_off = -1.1 #default = -1.1

score_cut_off = float(sys.argv[2])

min_atoms_in_sl = 4	#how many atoms superligands does at least have to have

print('superligand density: ', min_hits)
print('superligand cut off score: ', score_cut_off)

#*******************************************************************************


#---------------------------------------------
def gen_mol_list(cut_off):
	mol_list = []
	outdock = open('OUTDOCK', 'r')

	for i in outdock.readlines():
		if (i[0:6] == '     E') and (string.strip(i[72:79]) <> 'Total'):
			name = i[7:16]
			score = float( string.strip( i[71:79] ) )
			vdw = float( string.strip( i[37:48] ) )
			if score <> vdw:
				print('total score <> vdw score <-------------------------------------------')
				print(i)
				
			else:
				ratio = vdw/nhvy
				if ratio <= cut_off:
					mol_list.append(name)

					print(name, ratio, vdw, nhvy)

		else:
			try:
			nhvy = float( string.strip( i[39:41] ) ) #nhvy is in line before score
			except:
				continue
	outdock.close()
	return mol_list

#---------------------------------------------
def keep_pdb(mol_list,atom_cut_off):
	atom_dict = {}
	line_dict ={}
	found = False
	start = False
	poses = open('res.eel1', 'r')
	superlig = open('superligand.pdb', 'w')
	for i in poses.readlines():
		label = i[7:16]
		if i[0:6] == 'REMARK' and label in mol_list:
			#print 'consider', label
			label = i[7:16]
			#print label
			start = True
			found = True
		elif start and i.strip() <> 'TER' and i[13] <> 'H': # do not keep hydrogen atoms
			x = float(i[30:38])
			y = float(i[38:46])
			z = float(i[46:54])
			#print i
			#print x,y,z
			if len(atom_dict) == 0: #populate list
				atom_dict[x,y,z] = 1
				line_dict[x,y,z] = i
				found = True  #found at least one molecule that could be docked and meets cut off criteria
			else:
				#print atom_list
				atom_list = atom_dict.keys()
				for entry in atom_list: #compare new atom with every atom in atom list
					reject = False
					dsquare = math.pow(atom_cut_off,2)
					#print dsquare, 'here'
					distance = math.pow(entry[0]-x,2) + math.pow(entry[1]-y,2) + math.pow(entry[2]-z,2)
					#print distance
					#print distance, 'distance'
					if distance < dsquare: # compare square of distance, not distance.. Auri claims that this is computationally faster	
						#atom is close to an existing atom -> reject
						reject = True
						#mark how found close atom was found
						atom_dict[entry] = atom_dict[entry] + 1
						#print reject
						#print 'reject atom'
						break #no need to compare the atom any further
				if not reject:
					atom_dict[x,y,z] = 1
					line_dict[x,y,z] = i
					
		elif i.strip() == 'TER':
			start = False

	poses.close()

	sl_atom_counter = 0

	for key in line_dict.keys():
		if atom_dict[key] > min_hits:
			superlig.write(line_dict[key][0:55] + '\n')
			sl_atom_counter = sl_atom_counter + 1

	superlig.close()

	if not found or sl_atom_counter < min_atoms_in_sl :
		print('No atoms in superligand for this PDB Assuming pocket was too tight for docking. Substituting superligand with original ligand')
		command = 'cp ../xtal-lig.pdb superligand.pdb'
		os.system(command)

	#AK and AS run correctPDB.py after superligand, not sure why this is needed


#---------------------------------------------

#MAIN

mol_list = gen_mol_list(score_cut_off)
keep_pdb(mol_list, atom_cut_off)


#Spicoli has problems with F,Br and I -> change to C
super_lig_file = open("superligand.pdb", "r")
lines = super_lig_file.readlines()
super_lig_file.close()
super_lig_new = open("superligand.pdb", "w")


counter = 1000

for line in lines:
	line = 'ATOM   ' + str(counter) + '  CA  ASP  ' + str(counter)  + line[26:54] + '\n'
	counter = counter + 1

	super_lig_new.write(line)

super_lig_new.close()







		
 

