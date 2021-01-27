#!/usr/bin/env python3.6

# Cmd line version

import os, sys, math
import csv


from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFreeSASA, RDConfig, rdmolops, rdMolTransforms, rdGeometry, Descriptors3D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from argparse import ArgumentParser


from math import degrees
import time
import numpy

import openbabel

#--------------------------------------------------------------------------------------------------------------------


"""
radius dict taken from Tsai J, Taylor R, Chothia C, Gerstein M. J Mol Biol. 1999 Jul 2;290(1):253-66. https://www.sciencedirect.com/science/article/pii/S0022283699928292?via%3Dihub
RNA radii taken from NucProt: Gerstein, Voss J. Mol. Biol.(2005)346, 477–492, http://geometry.molmovdb.org/files/libproteingeometry/data/NucOr.atom-vols.dat http://geometry.molmovdb.org/files/libproteingeometry/data/NucProt.atom-defs.dat

O2H0b (1.62 Å) for ribose O3' and O5'-oxygen atoms, while O2H0s (1.50 Å) are assigned to O4'
"""
radius_dict_protor = {'C3H0': 1.61, 'C3H1':1.76,'C4H1':1.88, 'C4H2':1.88, 'C4H3':1.88, "N2H0": 1.64, 'N3H0':1.64, 'N3H1':1.64, 'N3H2':1.64,'N4H3':1.64, 'O1H0':1.42, 'O2H0': 1.62, 'O2H0s':1.50, 'O2H0b':1.62, 'O2H1':1.46, 'S2H0':1.77, 'S2H1':1.77,'C4H4':1.88, 'O2H0':1.46, 'P4H0': 1.82, 'P5H1':1.8, 'P3H0':1.8, "Mg2H2": 1.73, "Zn0H0":0.139, "F1H1": 1.0, "FE": 1.70, "MG": 1.59}
#radius dict from Jerry Tsai,RobinTaylor, Cyrus Chothia Mark Gerstein Mol. Biol.(1999)290, 253-266 (https://reader.elsevier.com/reader/sd/pii/S0022283699928292?token=5F39ED134A3E204BF8CDE9593BDFCE0C789DBB171F0FEF3ADFC8AC6B21C17D60A4EFAFAD625A7A06E4851D4FC3F07799)
radius_dict_bondi = {'C3H0': 1.61, 'C3H1':1.70,'C4H1':1.70, 'C4H2':1.70, 'C4H3':1.70, 'N3H0':1.55, 'N3H1':1.55, 'N3H2':1.55,'N4H3':1.55, 'O1H0':1.52,'O2H1':1.52, 'S2H0':1.80, 'S2H1':1.80,'C4H4':1.70, 'O2H0':1.52, 'P5H1':1.8, 'P3H0':1.8, "Mg2H2": 1.73, "Zn0H0":0.139}

#--------------------------------------------------------------------------------------


nuc_dict = {}
nuc_dict['A'] = ["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"]
nuc_dict['G'] = ["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"]
nuc_dict['C'] = ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
nuc_dict['U'] = ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
nuc_dict['T'] = ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C7"]

one_letter_prot ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q','ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y','ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

one_letter_nuc = ["A", "C", "G", "U", "T"]

#list of atom names in ribose
rib_list = ["O1P", "O2P", "P", "05'", "C5'", "C4'", "O4'", "C1'", "C2'", "C3'", "O2'", "O3'"]


#--------------------------------------------------------------------------------------

RNA_atom_types = {
"A" : {"N1":"N2H0" ,"C2":"C3H1", "N3":"N2H0", "C4":"C3H0", "C5": "C3H0", "C6":"C3H0", "N6": "N3H2", "N7":"N2H0", "C8":"C3H1", "N9":"N3H0"},
"G" : {"N1":"N3H1" ,"C2":"C3H0", "N2":"N3H2", "N3":"N2H0", "C4":"C3H0", "C5":"C3H0", "C6":"C3H0", "O6":"O1H0", "N7":"N2H0", "C8":"C3H1", "N9":"N3H0"},
"C" : {"N1":"N3H0", "C2":"C3H0", "O2":"O1H0", "N3":"N2H0", "C4":"C3H0", "N4":"N3H2", "C5":"C3H1", "C6":"C3H1"},
"U" : {"N1":"N3H0", "C2":"C3H0", "O2":"O1H0", "N3":"N2H0", "C4":"C3H0", "O4":"O1H0", "C5":"C3H1", "C6":"C3H1"},
"T" : {"N1":"N3H0", "C2":"C3H0", "O2":"O1H0", "N3":"N2H0", "C4":"C3H0", "O4":"O1H0", "C5":"C3H0", "C6":"C3H1", "C7":"C4H3"},
"R" : {"O1P":"O2H0", "OP1":"O2H0", "O2P":"O2H1", "OP2":"O2H1", "OP3":"O2H1", "P":"P4H0", "O5'":"O2H0b", "C5'":"C4H2", "C4'":"C4H1", "O4'":"O2H0s", "C1'":"C4H1", "C2'":"C4H1", "C3'":"C4H1", "O2'":"O2H1", "O3'":"O2H0b"},
}
#OP3 is terminal oxygen?


#--------------------------------------------------------------------------------------------------------------------
#spicoli is buggy, if CL,BR,or I present in superligand (F seems to be fine) surface of atoms of superligand complex outside binding site changes!!!!
#TO DO: set all superligand atoms to C
#introduce check, that radius of first atom in prot file and sl complex is the same
#added code to calculate aromatic SA and charged SA

docking_scripts = os.environ['DockingScripts']

#*******************************************************************************
# MAIN help
#**********

description = "Script to carry out descriptor calculations of DrugPred_RNA"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

parser.add_argument("-prot",  action="store", dest="prot", help="PDB ID")
parser.add_argument("-id", "--id", type=str, action="store", dest="id", help="ID code of prot",required=True)
parser.add_argument("--debug",  action="store_true", dest="debug",default=False, help="Print debugging output")
#parser.add_argument("-metal",  action="store", dest="metal",default=False, help="Metal in binding site", nargs='*')
#parser.add_argument("-lig",  action="store", type = str, dest="lig",default=False, help="Ligand in binding site")
#parser.add_argument("-s","--sl_cut_off_score", action="store", dest="sl_cut_off", default='-1.1', help="Score cut off for superligand, default = -1.1")
#parser.add_argument("--descriptors", action="store", dest="des_type", default='new', help="old or new type of descriptors (modify SA for N/O)",choices=['new', 'old', 'test'])

args = parser.parse_args()
args_dict = vars(args)
locals().update(args_dict) #generates local variables from key / value pairs

print("args;", args_dict)


#-------------------------------------------------
def TransferProp(mol, OBAtomPropertyQuery, propname, prop):
	"""
	Transfers properties from OpenBabel to RDKit via the SetBoolProp() function.
	Sets a custom BoolProp that is not previously available in the RDKit atom class.
	Indexing is adjusted, as RDKit starts with 0 and OB starts with 1

	"""
	prop = str(prop.lower())
	if "bool" in prop:
		if OBAtomPropertyQuery: #e.g. obatom.IsHbondAcceptor()
			mol.GetAtomWithIdx(obatom.GetIdx()-1).SetBoolProp(propname, True)
		else:
			mol.GetAtomWithIdx(obatom.GetIdx()-1).SetBoolProp(propname, False)
	elif "str" or "string" or "float" in prop:
		if OBAtomPropertyQuery: #e.g. obatom.IsHbondAcceptor()
			mol.GetAtomWithIdx(obatom.GetIdx()-1).SetProp(propname, str(OBAtomPropertyQuery))
		else:
			mol.GetAtomWithIdx(obatom.GetIdx()-1).SetProp(propname, "Not Applicable")
	elif "int" or "Int" in prop:
		if OBAtomPropertyQuery: #e.g. obatom.IsHbondAcceptor()
			mol.GetAtomWithIdx(obatom.GetIdx()-1).SetIntProp(propname, OBAtomPropertyQuery)
		else:
			print("Not Applicable")
	else:
		return print("Property undefined")


def block_apolar_N_atoms(mol, conf, cplx_mol, cplx_conf, polar_list): #places two C-atoms 1.5 Å above and below aromatic N-atoms and exocyclic N and O atoms adjacent to aromatic rings.
	"""
	Method:
	The unit vector is found from two neighbouring atoms. This unit vector is scaled to the appropiate distance of the two new dummy atoms from the N-atom (2.4 Å)
	The two atoms are placed at the N-atom coord + and - the scaled unit vector. The cplx is input here because the atom properties of the prot is well defined,
	but not for the cplx (e.g. GetIsAromatic()) Blame the sanitize-function of RDKit. 

	"""
	scalar = 1.70
	blocking_element = Chem.Atom(6) #which element is used for blocking (6 = carbon, 9 = Fluoride)
	global number_c_atoms

	mw = Chem.RWMol(mol)
	mw_cplx = Chem.RWMol(cplx_mol)
	mw_conf = mw.GetConformer()
	mw_cplx_conf = mw_cplx.GetConformer()
	#print(mw, mw_conf)
	block_delta_SASA_total = 0.0
	for polar_idx in polar_list:
		#print("polar idx", polar_idx)
		atom = mol.GetAtomWithIdx(polar_idx)
		nbr_list = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
		idx = 0
		#print(atom, atom.GetSymbol(), atom.GetIsAromatic(), atom.IsInRing(), atom.GetDegree())
		if atom.GetSymbol() == "N" and atom.IsInRing() and atom.GetIsAromatic() and (atom.GetDegree() == 2 or atom.GetDegree() == 3): # This atom is inside a heterocycle, and will be used with its two neighboring atoms to determine the unit vector 
			#print("Found tryptophan or similar")
			idx = atom.GetIdx()
			#print("idx", idx)
			#print("Found aromatic N")
			atom_3d_point = conf.GetAtomPosition(idx)
			#print("pos", list(atom_3d_point))
			n = 0
			vector_list = []
			while n < 2: #need two neighbors to find cross product -> unit vector
				for nbr in atom.GetNeighbors():
					nbridx = nbr.GetIdx()
					#print("nbridx", nbridx)
					nbr_3d_point = conf.GetAtomPosition(nbridx)
					vector = nbr_3d_point - atom_3d_point
					#print(vector)
					#print("vector", vector, list(vector))
					vector_list.append(vector)
					n += 1
			unit_vector = rdGeometry.Point3D.CrossProduct(vector_list[0], vector_list[1]) #the unit vector is perpendicular to the two vectors
			#print("unit_vector", unit_vector, list(unit_vector))
			length_unit_vector = rdGeometry.Point3D.Length(unit_vector) #this is usually ~1.72 Å for a N in an aromatic six-membered ring
			#print("length:", length_unit_vector, "Å")
			scaled = scalar / length_unit_vector #let's figure out how much the unit vector needs to be scaled
			scaled_vector =[i * scaled for i in list(unit_vector)] #let's scale the vector
			scaled_vector_3d = rdGeometry.Point3D(scaled_vector[0], scaled_vector[1], scaled_vector[2]) #make a Point3D object. Easier to handle in RDKit
			A1 = atom_3d_point + scaled_vector_3d # First atom coordinate, over ring
			A2 = atom_3d_point - scaled_vector_3d # Second atom coordinate, under ring
			#print("A1", list(A1), A1,"A2", list(A2), A2)
			A1_idx = mw.AddAtom(blocking_element) # first added atom to prot
			B1_idx = mw_cplx.AddAtom(blocking_element) # first added atom to cplx
			#print("A1, B1", A1_idx, B1_idx)
			C1 = mw_conf.SetAtomPosition(A1_idx, A1) # change coordinate for first added atom
			C2 = mw_cplx_conf.SetAtomPosition(B1_idx, A1) # change coordinate for first added atom
			A2_idx = mw.AddAtom(blocking_element)
			B2_idx = mw_cplx.AddAtom(blocking_element)
			D1 = mw_conf.SetAtomPosition(A2_idx, A2)
			D1 = mw_cplx_conf.SetAtomPosition(B2_idx, A2)
			#print("A1_idx", A1_idx, "A2_idx", A2_idx)
			#print("A2_idx", A2_idx)
			
			#calculate SASA
			mw_radii = prot_radii + [radius_dict_protor["C4H4"]] + [radius_dict_protor["C4H4"]] #add two carbon radii for the blocking atoms
			mw_cplx_radii = prot_radii + sl_radii + [radius_dict_protor["C4H4"]] + [radius_dict_protor["C4H4"]]
			rdFreeSASA.CalcSASA(mw, mw_radii, opts=opts)
			rdFreeSASA.CalcSASA(mw_cplx,  mw_cplx_radii, opts=opts)

			#calculate blocked SASA
			mw_SASA = mw.GetAtomWithIdx(polar_idx).GetProp("SASA")
			mw_cplx_SASA = mw_cplx.GetAtomWithIdx(polar_idx).GetProp("SASA")
			block_delta_SASA = float(mw_SASA) - float(mw_cplx_SASA)

			#calculate unblocked SASA ("Regular SASA")
			unblock_SASA = prot.GetAtomWithIdx(polar_idx).GetProp("SASA")
			unblock_cplx_SASA = cplx.GetAtomWithIdx(polar_idx).GetProp("SASA")
			reg_delta_SASA = float(unblock_SASA) - float(unblock_cplx_SASA)

			#calculate the difference between the blocked vs unblocked SASA. This surface area will be appended to the hydrophobic area and subtracted from the polar donor acceptor area.
			blocked_SA = reg_delta_SASA - block_delta_SASA 
			block_delta_SASA_total += blocked_SA # this is the blocked surface area

			if blocked_SA / reg_delta_SASA >= 0.5: #if more than half of the surface area is converted from polar to hydrophobic area, this is now considered as hydrophobic
				number_c_atoms += 1
			mw.RemoveAtom((mw.GetNumAtoms()-1)) #remove the last atom. Repeat three times so blocking atoms above and under are removed both in prot and cplx
			mw_cplx.RemoveAtom((mw.GetNumAtoms()-1))
			mw.RemoveAtom((mw.GetNumAtoms()-1))
			mw_cplx.RemoveAtom((mw.GetNumAtoms()-1))


		#print(atom.GetPropsAsDict())
		#print(atom.GetProp("atom_name"))
		elif atom.GetProp("atom_name") in ["O4", "O6", "O2", "N4", "N6", "N2"]: #This is an exocyclic nitrogen or oxygen atom. The nearest neighbor + its two nearest neigbor will be used to determine the unit vector, but it's still the exocyclic which will be blocked
			#print("Gotem", atom.GetSymbol(), atom.GetProp("atom_name"))
			vector_list = []
			O_atom_3d_point = conf.GetAtomPosition(atom.GetIdx())
			#if atom.GetSymbol() == "O":
			#	print("Oxygen atom x", atom.GetSymbol(), O_atom_3d_point.x, O_atom_3d_point.y, O_atom_3d_point.z)
			#if atom.GetSymbol() == "N":
			#	print("Nitrogen atom x", atom.GetSymbol(), O_atom_3d_point.x, O_atom_3d_point.y, O_atom_3d_point.z)
			for nbratom in atom.GetNeighbors(): #1st neighbor, should be carbon in the aromatic ring
				nbridx = nbratom.GetIdx()
				nbr_3d_point_1 = conf.GetAtomPosition(nbridx)
				#print("1st nbr x", nbr_3d_point_1.x)
				#print(mol.GetAtomWithIdx(nbridx))
				nbr_atom = mol.GetAtomWithIdx(nbridx)
				for nbr in nbr_atom.GetNeighbors(): #2nd neigbors
					nbridx = nbr.GetIdx()
					#print("nbridx", nbridx)
					nbr_3d_point_2 = conf.GetAtomPosition(nbridx)
					#print(list(nbr_3d_point_2))
					if list(nbr_3d_point_2) != list(O_atom_3d_point): #comparing x-coordinate are not sufficient. Example: 3f4e, where C5 and and N6 in A87 have same x-coordinate
						#print("nbr atom x", nbr_3d_point_2.x)
						vector = nbr_3d_point_1 - nbr_3d_point_2
						#print(vector)
						#print("vector", vector, list(vector))
						vector_list.append(vector)

			unit_vector = rdGeometry.Point3D.CrossProduct(vector_list[0], vector_list[1]) #the unit vector is perpendicular to the two vectors
			#print("unit_vector", unit_vector, list(unit_vector))
			length_unit_vector = rdGeometry.Point3D.Length(unit_vector) #this is usually ~1.72 Å for a N in an aromatic six-membered ring
			#print("length:", length_unit_vector, "Å")
			scaled = scalar / length_unit_vector #let's figure out how much the unit vector needs to be scaled
			scaled_vector =[i * scaled for i in list(unit_vector)] #let's scale the vector
			scaled_vector_3d = rdGeometry.Point3D(scaled_vector[0], scaled_vector[1], scaled_vector[2]) #make a Point3D object. Easier to handle in RDKit
			A1 = O_atom_3d_point + scaled_vector_3d # First atom coordinate, over ring
			A2 = O_atom_3d_point - scaled_vector_3d # Second atom coordinate, under ring
			#print("A1", list(A1), A1,"A2", list(A2), A2)
			A1_idx = mw.AddAtom(blocking_element) # first added atom to prot
			B1_idx = mw_cplx.AddAtom(blocking_element) # first added atom to cplx
			#print("A1, B1", A1_idx, B1_idx)
			C1 = mw_conf.SetAtomPosition(A1_idx, A1) # change coordinate for first added atom
			C2 = mw_cplx_conf.SetAtomPosition(B1_idx, A1) # change coordinate for first added atom

			A2_idx = mw.AddAtom(blocking_element)
			B2_idx = mw_cplx.AddAtom(blocking_element)

			D1 = mw_conf.SetAtomPosition(A2_idx, A2)
			D1 = mw_cplx_conf.SetAtomPosition(B2_idx, A2)

			#calculate SASA
			mw_radii = prot_radii + [radius_dict_protor["C4H4"]] + [radius_dict_protor["C4H4"]] #add two carbon radii for the blocking atoms
			mw_cplx_radii = prot_radii + sl_radii + [radius_dict_protor["C4H4"]] + [radius_dict_protor["C4H4"]]
			rdFreeSASA.CalcSASA(mw, mw_radii, opts=opts)
			rdFreeSASA.CalcSASA(mw_cplx,  mw_cplx_radii, opts=opts)

			#calculate blocked SASA
			mw_SASA = mw.GetAtomWithIdx(polar_idx).GetProp("SASA")
			mw_cplx_SASA = mw_cplx.GetAtomWithIdx(polar_idx).GetProp("SASA")
			block_delta_SASA = float(mw_SASA) - float(mw_cplx_SASA)

			#calculate unblocked SASA ("Regular SASA")
			unblock_SASA = prot.GetAtomWithIdx(polar_idx).GetProp("SASA")
			unblock_cplx_SASA = cplx.GetAtomWithIdx(polar_idx).GetProp("SASA")
			reg_delta_SASA = float(unblock_SASA) - float(unblock_cplx_SASA)

			#calculate the difference between the blocked vs unblocked SASA. This surface area will be appended to the hydrophobic area and subtracted from the polar donor acceptor area.
			blocked_SA = reg_delta_SASA - block_delta_SASA 
			block_delta_SASA_total += blocked_SA # this is the blocked surface area
			if blocked_SA / reg_delta_SASA >= 0.5: #if more than half of the surface area is converted from polar to hydrophobic area, this atom is now considered as hydrophobic
				number_c_atoms += 1

			#the two latest added atoms (the blocking atoms) are removed from the structure, so that the blocked atoms surface area does not infer with later calculations
			mw.RemoveAtom((mw.GetNumAtoms()-1))
			mw_cplx.RemoveAtom((mw.GetNumAtoms()-1))
			mw.RemoveAtom((mw.GetNumAtoms()-1))
			mw_cplx.RemoveAtom((mw.GetNumAtoms()-1))
	if mw:
		#print("Num atoms:", mw.GetNumAtoms(), mol.GetNumAtoms())
		return mw, mw_conf, mw_cplx, mw_cplx_conf, block_delta_SASA_total

#-------------------------------------------------

def create_bs_seq(mol_in, bs_atom_list):
	binding_site_seq = []
	checked_res = []
	for idx in bs_atom_list:
		original_atom = mol_in.GetAtomWithIdx(idx)
		res_info = original_atom.GetPDBResidueInfo()
		original_res_num = res_info.GetResidueNumber()
		original_res_name = res_info.GetResidueName()
		if original_res_num not in checked_res:
			checked_res.append(original_res_num)
			if original_res_name.strip() in one_letter_nuc:
				one_letter = original_res_name.strip()
			elif original_res_name.strip() in one_letter_prot:
				one_letter =one_letter_prot[original_res_name.strip()]
			binding_site_seq.append(one_letter)

	return binding_site_seq

#----------------------------------------------------

def create_bs_mol(mol_in, mol_in_conf, bs_atom_list):
	"""
	Creates a mol of the binding site atoms, e.g. the atoms where deltaSASA is larger than zero.
	Validates f_C_atoms
	"""
	empty_mol = Chem.Mol() #create empty mol
	mol_rw = Chem.RWMol(empty_mol) #make empty mol editable. This will be our output
	conformer =  Chem.Conformer(len(bs_atom_list)) # create conformer for new mol. Otherwise, unable to add xyz
	for idx in bs_atom_list:
		#get original atom info
		original_atom = mol_in.GetAtomWithIdx(idx)
		res_info = original_atom.GetPDBResidueInfo()
		original_element =  original_atom.GetAtomicNum()
		original_res_num = res_info.GetResidueNumber()
		original_res_name = res_info.GetResidueName()
		original_xyz =  list(mol_in_conf.GetAtomPosition(idx)) #
		#set new atom, info
		new_xyz = rdGeometry.Point3D(original_xyz[0], original_xyz[1], original_xyz[2]) #create Point3D object which can be added to conformer
		new_atom_idx = mol_rw.AddAtom(Chem.Atom(original_element)) #add atom, returns new atom idx
		new_atom = mol_rw.GetAtomWithIdx(new_atom_idx)
		#print(new_atom_idx)
		new_res_inf = Chem.AtomPDBResidueInfo(res_info.GetName(), serialNumber=res_info.GetSerialNumber(),
                                       altLoc=res_info.GetAltLoc(), residueName=res_info.GetResidueName(),
                                       residueNumber=res_info.GetResidueNumber(), chainId=res_info.GetChainId(),
                                       insertionCode=res_info.GetInsertionCode(), occupancy=res_info.GetOccupancy(), \
                                       tempFactor=res_info.GetTempFactor(), isHeteroAtom=res_info.GetIsHeteroAtom(),
                                       secondaryStructure=res_info.GetSecondaryStructure(),
                                       segmentNumber=res_info.GetSegmentNumber())
		new_atom.SetMonomerInfo(new_res_inf) #invoking this line creates trouble
		new_atom = mol_rw.GetAtomWithIdx(new_atom_idx)
		conformer.SetAtomPosition(new_atom_idx, new_xyz) #assign new xyz to new mol
	


	mol_out = mol_rw.GetMol()
	mol_out_conf = mol_out.AddConformer(conformer, assignId = True) #merges conformer with mol

	#mi  =  Chem.AtomPDBResidueInfo()
	#mi.SetResidueName('MOL')
	#mi.SetResidueNumber(1)
	#mi.SetOccupancy(0.0)
	#mi.SetTempFactor(0.0)
	#[a.SetMonomerInfo(mi)  for  a  in  mol_out.GetAtoms()]
	

	return mol_out, mol_out_conf


#-------------------------------------------------


def determine_radii_sl(molecule): #this function should only be applied for the superligand atoms
	radii = []
	for atom in molecule.GetAtoms():
		radii.append(radius_dict_protor["C4H4"])
	return radii

def determine_radii(molecule):
	radii=[]
	atom_type_list = []
	pdb_name_list = []
	#molecule.UpdatePropertyCache(strict=True) #Recalculates valence states of atoms. See issue #1596, https://github.com/rdkit/rdkit/issues/1596. If used, KDOPS and IDH gives errors on valency
	for atom in molecule.GetAtoms():
		info = atom.GetPDBResidueInfo()
		residue = info.GetResidueName().strip()

		atom_name = atom.GetProp("atom_name")
		#atom_type = atom.GetSymbol() + str(atom.GetProp("HV")) + 'H' + str(atom.GetProp("IHC"))
		prot_atom_type = atom.GetSymbol() + str(atom.GetTotalDegree()) + "H" + str(atom.GetNumImplicitHs())
		#print(prot_atom_type)
		#print(atom_type)
		if atom_name.strip() in RNA_atom_types["R"]:
			atom_type_list.append(atom_name)
			pdb_name_list.append(RNA_atom_types["R"][atom_name])
			radii.append(radius_dict_protor[RNA_atom_types["R"][atom_name]])
			atom.SetDoubleProp("radius", radius_dict_protor[RNA_atom_types["R"][atom_name]])
			atom.SetProp("protor_name", RNA_atom_types["R"][atom_name])
		elif residue in RNA_atom_types:
			#print(residue, atom_name)
			if atom_name in RNA_atom_types[residue]:
				#print(RNA_atom_types[residue][atom_name])
				atom_type_list.append(atom_name)
				pdb_name_list.append(RNA_atom_types[residue][atom_name])
				radii.append(radius_dict_protor[RNA_atom_types[residue][atom_name]])
				atom.SetDoubleProp("radius", radius_dict_protor[RNA_atom_types[residue][atom_name]])
				atom.SetProp("protor_name", RNA_atom_types[residue][atom_name])
		elif residue in mod_res_dic:

			mimick =  mod_res_dic[residue]
			if mimick in RNA_atom_types:
				#print("MODRES found:", residue, atom_name, "MIMICKS:", mod_res_dic[residue])
				if atom_name in RNA_atom_types[mimick]:
					#print(RNA_atom_types[mimick][atom_name])
					atom_type_list.append(atom_name)
					pdb_name_list.append(RNA_atom_types[mimick][atom_name])
					radii.append(radius_dict_protor[RNA_atom_types[mimick][atom_name]])
					atom.SetDoubleProp("radius", radius_dict_protor[RNA_atom_types[mimick][atom_name]])
					atom.SetProp("protor_name", RNA_atom_types[mimick][atom_name])
				else:
					#print("Atom not assigned:", atom_name)
					atom.SetDoubleProp("radius", 0.0 )
					atom.SetProp("protor_name", "Unassigned")
					radii.append(Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()))
			else:
				#print("Atom not assigned:", atom_name)
				atom.SetDoubleProp("radius", 0.0 )
				atom.SetProp("protor_name", "Unassigned")
				radii.append(Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()))
		elif prot_atom_type in radius_dict_protor:
			radii.append(radius_dict_protor[prot_atom_type])
			atom.SetProp("protor_name", prot_atom_type)
			atom.SetDoubleProp("radius", radius_dict_protor[prot_atom_type])
		elif atom_name in radius_dict_protor:
			radii.append(radius_dict_protor[atom_name])
			atom.SetProp("protor_name", atom_name)
			atom.SetDoubleProp("radius", radius_dict_protor[atom_name])
		else:
			#print (atom_type, 'not found', atom.GetIdx(), Chem.GetPeriodicTable().GetRcovalent(atom.GetAtomicNum()))
			#use VDW radius
			radii.append(Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()))
			#print("Default radius used:", atom_name, residue, Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()))
			atom.SetProp("protor_name", atom_name)
			atom.SetDoubleProp("radius", Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()))
		#print("atom_type_debug", atom_type, "idx", atom.GetIdx())
	return radii, atom_type_list, pdb_name_list

def get_atom_type(prot_cofac): # iterates over the _cofac.pdb file to retrieve the atom type string (eg. CG, N, NZ, OE1, OE2). Neither RDKit nor OpenBabel adheres to the PDB standard when retrieving these names.
	pdb_dic = {}
	idx_x_dic = {}
	path = "./" + prot_cofac
	cofac_file = open(path, "r")
	cofac_lines = cofac_file.readlines()
	pdb_idx = 0
	for cofac_line in cofac_lines:
		pdb_idx += 1
		atom_name = cofac_line[12:16].strip() #atom type
		#print("atom_name", atom_name, pdb_idx)
		pdb_dic[str(pdb_idx)] = atom_name
		idx_x_dic[str(cofac_line[30:38].strip())] = pdb_idx
	return pdb_dic, idx_x_dic

def attribute_atom_type(pdb_dic, mol):
	for atom in mol.GetAtoms():
		idx = atom.GetIdx() + 1
		if idx <= prot.GetNumAtoms():
			atom.SetProp("atom_name", pdb_dic[str(idx)])
		else:
			atom.SetProp("atom_name", "Undefined")

#-------------------------------------------------

### Version

print("RDKit version:", rdBase.rdkitVersion)


### Declare list of hydrophobic residues and hydrophobicity index dictionary

hydrophobic_reslist = ['ALA', 'GLY', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'PRO']
hydrophobicity_indices = {'ALA':1.8, 'ARG':-4.5, 'ASN':-3.5, 'ASP':-3.5, 'CYS':2.5, 'GLN':-3.5, 'GLU':-3.5, 'GLY':-0.4, 'HIS':-3.2, 'ILE':4.5, 'LEU':3.8, 'LYS':-3.9, 'MET':1.9, 'PHE':2.8, 'PRO': -1.6, 'SER':-0.8, 'THR':-0.7, 'TRP':-0.9, 'TYR':-1.3, 'VAL':4.2}
aa_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
charged_atoms = ['OE1', 'OE2', 'NZ','OD1', 'OD2','NE','NH1','NH2','ZN','MG','CA','MN','FE','NI','OP1','OP2'] 

### Declare values for grid spacing and probe radius during surface area calculations


#get values from db
global debug

print(prot)

mod_res_dic ={}

code_list = []
with open("../../../input_csvs/modified_bases_rna_dna.csv", newline = '\n') as csvfile:
	mod_res_dic_t =  [row for row in csv.reader(csvfile, delimiter= "\t")]

for row in mod_res_dic_t:
	mod_res_dic[row[0]] = row[1].upper()



#print(mod_res_dic)


use_cleaned_input = False

try:
	prot_cofac = prot + '_cofac.pdb'
	prot = Chem.MolFromPDBFile(prot_cofac, sanitize = False)

	#prot_san = Chem.MolFromPDBFile(prot_cofac)
	print(prot)
	print ('end)')
	protconf = prot.GetConformer()

except: 
	print ("Protein file unreadable")
	sys.exit(0)



### Call in protein-pdb via OpenBabel, get properties

aromatic_idx_list = []
acc_list = []
don_list = []
hyb_dic = {}
type_dic = {}

command = 'cat ' + prot_cofac + ' ./testing/superligand.pdb > sl_' + prot_cofac
print (command)
os.system(command)

pdb_lines, idx_x_dic = get_atom_type(prot_cofac)


cplx = Chem.MolFromPDBFile('sl_' + prot_cofac,sanitize=False)
cplxconf = cplx.GetConformer()

print("cplx succesfully read")
#### Call in protein


#print("pdb_lines", pdb_lines)

attribute_atom_type(pdb_lines, prot)

obConversion = openbabel.OBConversion()
obConversion.SetInFormat("pdb") # sets the format for the read file
OBmol = openbabel.OBMol() #makes an openbabel molecule

babelprot = obConversion.ReadFile(OBmol, prot_cofac) #reads in pdb-file
OBmol.UnsetAromaticPerceived() # unflags all the aromatic flags, so that OpenBabel properly recalculates the atoms aromatic property

failure = 0
#print(OBmol.NumConformers())
for obatom in openbabel.OBMolAtomIter(OBmol):
	failure += 1#
	if obatom.GetX() != protconf.GetAtomPosition(obatom.GetIdx()-1).x: #Compares the x-coordinates
		print ("idx:",obatom.GetIdx(), "x-coord", obatom.GetX(), "atom mismatch", protconf.GetAtomPosition(obatom.GetIdx()-1).x)
		print("index", obatom.GetIdx(), idx_x_dic[str(obatom.GetX())])
		#print(idx_x_dic[str(obatom.GetX())], protconf.GetAtomPosition(idx_x_dic[obatom.GetX()]).x)
		#print("obidx:", obatom.GetIdx(), "x-coords:", obatom.GetX(), protconf.GetAtomPosition(obatom.GetIdx()-1).x)
		sys.exit()
	else:
		TransferProp(prot, obatom.IsHbondDonor(), "HBondDonor", "Bool")
		TransferProp(prot, obatom.IsHbondAcceptor(), "HBondAcceptor", "Bool")
		TransferProp(prot, obatom.GetPartialCharge(), "Partial charge", "Float")
		TransferProp(prot, obatom.GetFormalCharge(), "Formal charge", "Float")
		TransferProp(prot, obatom.GetResidue().GetName(), "residue_name", "String")
		TransferProp(prot, obatom.GetResidue().GetNum(), "Residue number", "String")
		TransferProp(prot, obatom.GetType(), "Atom type", "String")
		TransferProp(prot, obatom.GetValence(), "Valence", "Int")
		TransferProp(prot, obatom.ImplicitHydrogenCount(), "IHC", "Int")
		TransferProp(prot, obatom.ExplicitHydrogenCount(), "EHC", "Int")
		TransferProp(prot, obatom.GetHvyValence(), "HV", "Int")
		TransferProp(prot, obatom.GetImplicitValence(), "IV", "Int")#
		TransferProp(cplx, obatom.IsHbondDonor(), "HBondDonor", "Bool")
		TransferProp(cplx, obatom.IsHbondAcceptor(), "HBondAcceptor", "Bool")
		TransferProp(cplx, obatom.GetPartialCharge(), "Partial charge", "Float")
		TransferProp(cplx, obatom.GetFormalCharge(), "Formal charge", "Float")
		TransferProp(cplx, obatom.GetResidue().GetName(), "residue_name", "String")
		TransferProp(cplx, obatom.GetResidue().GetNum(), "Residue number", "String")
		TransferProp(cplx, obatom.GetType(), "Atom type", "String")
		TransferProp(cplx, obatom.GetValence(), "Valence", "Int")
		TransferProp(cplx, obatom.ExplicitHydrogenCount(), "EHC", "Int")
		TransferProp(cplx, obatom.GetHvyValence(), "HV", "Int")
		TransferProp(cplx, obatom.GetImplicitValence(), "IV", "Int")#
		#obatom_formal_charge = obatom.GetFormalCharge()
		#obatom_partial_charge =obatom.GetPartialCharge()
		#print(obatom_formal_charge, obatom_partial_charge)#
		if obatom.IsAromatic():
			prot.GetAtomWithIdx(obatom.GetIdx()-1).SetIsAromatic(True)
		else:
			prot.GetAtomWithIdx(obatom.GetIdx()-1).SetIsAromatic(False)


	

#read superligand
prot_radii, prot_type_list, pdb_name_list = determine_radii(prot)

#print("len_prot_radii", len(prot_radii))
#print("number of atom", prot.GetNumAtoms())

if len(prot_radii) != prot.GetNumAtoms():
	print("Incorrect assignment of radii")


#print(prot_radii)

superligand = Chem.MolFromPDBFile('./testing/superligand.pdb',proximityBonding=False,sanitize=False)

superligandconf = superligand.GetConformer()

#position = superligand.GetConformer().GetAtomPosition(5)


if not superligand: #superligand is None if reading failed
	print("Superligand file undefined")
	sys.exit(0)

sl_radii = determine_radii_sl(superligand)

#print("len sl_radii", len(sl_radii))
#print("num_sl_atoms", superligand.GetNumAtoms())

print("##############")


vol = AllChem.ComputeMolVolume(superligand) # core dump error


#### Call in protein-superligand complex

cplx_radii = prot_radii + sl_radii
cplx_block_radii = prot_radii + sl_radii
attribute_atom_type(pdb_lines, cplx)


print("Num atoms",cplxconf.GetNumAtoms())

#### generate surfaces 

#method = rdFreeSASA.SASAAlgorithm.vol
opts = rdFreeSASA.SASAOpts()
opts.algorithm = rdFreeSASA.LeeRichards
opts.probeRadius = 1.0

surf_prot = rdFreeSASA.CalcSASA(prot, prot_radii,opts=opts)

print("surf_prot:", surf_prot)

cplx_surf = rdFreeSASA.CalcSASA(cplx, cplx_radii,opts=opts)
print("cplx_surf:", cplx_surf)

sl_surf = rdFreeSASA.CalcSASA(superligand, sl_radii,opts=opts)
print("superligand surface:", sl_surf)
	
#print("BLOCK SURF", surf_prot_block, surf_cplx_block, surf_prot_block - surf_cplx_block)
print("UNBLOCK SURF", surf_prot, cplx_surf, surf_prot - cplx_surf)

### Calculate any SASA changes for every atom in protein when converted to complex

backbone_atoms = ['N','CA','C','O']
num_apolar_res = 0.0 # Float required, because otherwise division removes decimal points
num_res = 0.0
hydrophobic_sasa_total = 0.0
hydrophilic_sasa_total = 0.0
sasa_total = 0.0
aromatic_sasa_total = 0.0
charged_sasa_total = 0.0
sum_hydrophobicity_index = 0.0
don_acc_sasa = 0.0
aromatic_N_O_sasa = 0.0
c_aliphatic_sasa = 0.0
all_aliphatic_sasa = 0.0
number_c_atoms = 0.0
no_sl_atoms = 0.0



haa_covered_res = []
hiaa_covered_res = []

#----------------------------------------------------------------------- RDKit SASA Change start

prot_surface_area = 0.0


total_delta_sasa = 0.0
number_c_atoms = 0.0
hydrophobic_sasa_total = 0.0
c_aliphatic_sasa = 0.0
all_aliphatic_sasa = 0.0
aromatic_sasa_total = 0.0
sasa_total = 0.0
bs_atom = 0.0
dsasa_3 = 0.0 #Sum of all superligand SASA values in complex
buried_sl_atom = 0.0 #old name: exp_sl_atoms
surface_atoms = 0.0
buried_sl_atoms = 0.0

all_res = []
all_res_name = []
all_res_num = []
all_idx = []
all_chain = []
haa_covered_res = []
hiaa_covered_res = []
#binding_site_seq = []
index_list = []
polar_residue_list = []
bs_atom_list = []
bs_res_list = []
atom_idx_by_res = {}

cplx_surface_area = 0.0
hydrophobic_count = 0

monomer_info  =  Chem.AtomPDBResidueInfo()

for atom in prot.GetAtoms():
	info = atom.GetPDBResidueInfo()
	res_name = info.GetResidueName().strip()
	res_num = str(info.GetResidueNumber())
	all_res.append(res_num)
	all_res_name.append(info.GetResidueName())
	all_res_num.append(info.GetResidueNumber())
	#all_chain.append(info.GetChainId())
	#print(info.GetResidueNumber())
	all_idx.append(atom.GetIdx())
	protatom_SASA = atom.GetProp("SASA") #use this if you want the unblocked N-atom SASA
	#print(atom.GetIdx() + 1, atom.GetProp("residue_name"), atom.GetProp("atom_name"), atom.GetProp("protor_name"), atom.GetProp("radius"))
	#if atom.GetSymbol() == "O":
		#print(protconf.GetAtomPosition(atom.GetIdx()).x, protatom_SASA)
	#protatom_SASA =  prot_block.GetAtomWithIdx(atom.GetIdx()).GetProp("SASA") # retrieves the SASA info from the blocked mol. This mol doesn't play nice with the TransferProp-function
	
	if protconf.GetAtomPosition(atom.GetIdx()).x == cplxconf.GetAtomPosition(atom.GetIdx()).x:
		pass
	else:
		print("Atom mismatch")
	if res_name in one_letter_prot:
		res_name = one_letter_prot[res_name]
	atom_type = atom.GetProp("atom_name")
	prot_idx = atom.GetIdx()

	#print(prot_idx+1, atom.GetIsAromatic(), atom.GetProp("Valence"), atom.GetProp("IHC"), atom.GetProp("HV"), atom.GetProp("IV"))


	prot_surface_area = prot_surface_area + float(protatom_SASA)
	cplxatom_SASA = cplx.GetAtomWithIdx(prot_idx).GetProp("SASA") #use this one for the unblocked SASA
	#cplxatom_SASA= cplx_block.GetAtomWithIdx(atom.GetIdx()).GetProp("SASA") #use this one for the blocked SASA
	resid = res_name + res_num
	monomer_info.SetTempFactor(-20)
	#atom_name = '{0:>2s}{1:<2d}'.format(atom.GetSymbol(), 1)
	#monomer_info.SetName(atom_name)

	if protatom_SASA != cplxatom_SASA and (res_name in one_letter_nuc or one_letter_prot) :
		bs_atom += 1
		bs_atom_list.append(prot_idx)
		#print(res_name)
		bs_res_idx = info.GetResidueNumber()
		bs_res_list.append(bs_res_idx)
		#print(atom.GetIdx())
		binding_site_res = res_name
		prot.SetBoolProp("SASA_Changed", True)
		delta_SASA = float(protatom_SASA) - float(cplxatom_SASA)
		#print("DELTA SASA", delta_SASA)
		#print("DELTASASA1", delta_SASA, "IDX", atom.GetIdx())
		total_delta_sasa = total_delta_sasa + delta_SASA
		#monomer_info.SetTempFactor(delta_SASA)
		#atom.SetMonomerInfo(monomer_info)
		#print("Delta_SASA:", atom.GetIdx()+1, atom.GetSymbol(), str(delta_SASA)[0:8])
		if not atom.GetIsAromatic():
			all_aliphatic_sasa = all_aliphatic_sasa + delta_SASA
		if atom.GetSymbol() in ["C", "P", "S"]: #Carbon atom, hydrophobic
			number_c_atoms = number_c_atoms + 1 
			#print("Carbon atom 1")
			hydrophobic_sasa_total = hydrophobic_sasa_total + delta_SASA
			if not atom.GetIsAromatic():
				c_aliphatic_sasa = c_aliphatic_sasa + delta_SASA


		if atom.GetSymbol() == "N" and atom.IsInRing() and atom.GetIsAromatic() and (atom.IsInRingSize(6) or atom.IsInRingSize(5)): #ring N (bases + His + Trp, but not PRO)
			if atom.GetImplicitValence() == 0 and atom.GetDegree() == 3 and len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "C"]) == 3: #three C neighbours -> hydrophobic
				#print("neigh list", [x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "C"])
				#print("center atom idx", atom.GetIdx())
				hydrophobic_sasa_total = hydrophobic_sasa_total + delta_SASA
				number_c_atoms = number_c_atoms + 1 #is basically an C atom
				#print("Carbon atom 2")
			else:
				polar_residue_list.append(atom.GetIdx())
				don_acc_sasa = don_acc_sasa + delta_SASA
				# This is where the polar_segment function was
				#vectorcalc


		elif (atom.GetSymbol() == "N" or atom.GetSymbol() == "O") :  #N or O
			#print(atom.GetSymbol())
			don_acc_sasa = don_acc_sasa + delta_SASA
			#print("N or O")
			for nbr in atom.GetNeighbors():#get neighbour, needed for bases
				break #only need one atom
			if (nbr.IsInRing() and (nbr.IsInRingSize(6) or nbr.IsInRingSize(5) or nbr.IsInRingSize(4))) and atom.GetImplicitValence() != 1:#we are dealing with a base => -NH2 or =O
				polar_residue_list.append(atom.GetIdx())

			elif nbr.IsInRing() and (nbr.IsInRingSize(6) or nbr.IsInRingSize(5) or nbr.IsInRingSize(4)) and atom.GetSymbol() == "O":
				polar_residue_list.append(atom.GetIdx())

		elif atom.GetSymbol() == "P": #Phosphor atoms are usually inside phosphates and shielded from any interactions with the superligand. Therefore, they are basically c-atoms ("PATTY")
			hydrophobic_sasa_total += delta_SASA
			number_c_atoms += 1
		elif atom.GetAtomicNum() > 8: #all elements over oxygen added to polar surface area, e.g. sulphur, metals etc.
			don_acc_sasa = don_acc_sasa + delta_SASA 


		if atom.GetIsAromatic():
			aromatic_sasa_total = aromatic_sasa_total + delta_SASA
		#print("atom type:", atom.GetProp("Atom type").strip(), type(atom.GetProp("Atom type").strip()) )
		if atom_type in charged_atoms and atom.GetSymbol() !=  "C" : #CA only for metal, not Carbon 
			#print("Caught charged atom")
			charged_sasa_total = charged_sasa_total + delta_SASA
		if resid not in all_res:
		#	all_res.append(resid)
			num_res += 1
		if res_name in hydrophobic_reslist and atom_type not in backbone_atoms:
			if resid not in haa_covered_res:
				#print("caught hydrophobic res")
				num_apolar_res += 1
				haa_covered_res.append(resid)
		if res_name in aa_residues and atom_type not in backbone_atoms:
			if resid not in hiaa_covered_res:
				sum_hydrophobicity_index = sum_hydrophobicity_index + hydrophobicity_indices[resname]
				hiaa_covered_res.append(resid)


binding_site_seq = create_bs_seq(prot, bs_atom_list)
print(binding_site_seq)

# Create a binding site mol

bs_mol, bs_mol_conf = create_bs_mol(prot, protconf, bs_atom_list)
print("Size:", bs_mol.GetNumAtoms())
Chem.rdmolfiles.MolToPDBFile(bs_mol, id + "_" + str(opts.probeRadius) + "A_bs.pdb")

# Create a binding site mol composed of residues which contain atoms in direct contact with sl

res_idx_tuple = list(zip(all_res, all_idx, all_res_name, all_res_num))

all_res_mol_idx = []
flagged_res_name = []
flagged_res_num = []

for res, bs_idx, res_name, num in res_idx_tuple:
	#print(res, bs_idx)
	if str(res) in bs_res_list:
		#print("Atom number " + str(bs_idx) +" is in a residue")
		all_res_mol_idx.append(bs_idx)
		flagged_res_name.append(res_name)
		flagged_res_num.append(num)
print(all_res_mol_idx)

bs_mol_res, bs_mol_res_conf = create_bs_mol(prot, protconf, all_res_mol_idx)
print("Size:", bs_mol_res.GetNumAtoms())
Chem.rdmolfiles.MolToPDBFile(bs_mol_res, id + "_" + str(opts.probeRadius) + "A_bs_res.pdb")
block = Chem.rdmolfiles.MolToPDBBlock(bs_mol_res)

print(len(flagged_res_name))
print(len(flagged_res_num))




print("Block START")

#print("don_acc_sasa", don_acc_sasa)
#print("hsa_t", hydrophobic_sasa_total)
prot_block, protconf_block, cplx_block, cplxconf_block, delta_SASA_block_unblock = block_apolar_N_atoms(prot, protconf, cplx, cplxconf, polar_residue_list)

print("pre-block don_acc_sasa", don_acc_sasa, "pre-block hsa_t", hydrophobic_sasa_total, "pre-block don_acc_sasa",  don_acc_sasa )
hydrophobic_sasa_total += delta_SASA_block_unblock
don_acc_sasa -= delta_SASA_block_unblock
print("post-block don_acc_sasa", don_acc_sasa, "post-block hsa_t", hydrophobic_sasa_total,  "pre-block don_acc_sasa",  don_acc_sasa )

print("Block END")


print("Binding site seq:", str("".join(binding_site_seq)), "index_list:", " ".join(index_list), "Length:", "Seq length:", len(binding_site_seq))
print("Polar residues", polar_residue_list, "length", len(polar_residue_list))

bs_seq = str("".join(binding_site_seq))
print(bs_seq)

pdb_file_title = str(id) +"_blocked_N_atoms.pdb"
print(pdb_file_title)


print("apolar_res:", num_apolar_res, "sum_hydrophobicity_index", sum_hydrophobicity_index)
print("total delta SASA:", (total_delta_sasa))
print("PROT - CPLX AREA", prot_surface_area - cplx_surface_area)
print("HYDROPHOBIC_SASA_TOTAL:", hydrophobic_sasa_total)



# CALCULATIONS FOR THE ENCLOSURE DESCRIPTOR
# Calculate SASA change for superligand

superligand_delta_sasa_total = 0.0 # This variable will give the total change in SASA of superligand when it is in complex with protein.
superligand_sasa_total = 0.0 # This variable will contain the total SASA of the superligand when not complexed with protein (can not accessed directly because of box)


#get Descriptors3D descriptors for the superligand

Descriptors3D_dic = {'PMI1' : Descriptors3D.PMI1, 'PMI2' : Descriptors3D.PMI2, 'PMI3' : Descriptors3D.PMI3, 'NPR1' : Descriptors3D.NPR1, 'NPR2' : Descriptors3D.NPR2, 'RadiusOfGyration' : Descriptors3D.RadiusOfGyration, 'InertialShapeFactor' :  Descriptors3D.InertialShapeFactor, 'Eccentricity' : Descriptors3D.Eccentricity, 'Asphericity' : Descriptors3D.Asphericity, 'SpherocityIndex' : Descriptors3D.SpherocityIndex}
Descriptors3D_values= {}

for prop in Descriptors3D_dic:
	function = Descriptors3D_dic[prop]
	x = function(superligand)
	Descriptors3D_values[prop] = x


prot_num = prot.GetNumAtoms()

for supligatom in superligand.GetAtoms():
	sl_idx = supligatom.GetIdx()
	sl_atom = superligand.GetAtomWithIdx(sl_idx)
	orig_sl_info = sl_atom.GetPDBResidueInfo()
	no_sl_atoms += 1
	sl_coord = superligandconf.GetAtomPosition(sl_idx).x
	#cplx_sl_coord = cplxconf.GetAtomPosition(prot_num + sl_idx).x
	slatom_SASA = supligatom.GetProp("SASA")
	if slatom_SASA != 0:
		surface_atoms += 1
	else:
		buried_sl_atoms += 1

	cplx_sl_SASA = cplx.GetAtomWithIdx(prot_num + sl_idx).GetProp("SASA")
	dsasa_3 += float(cplx_sl_SASA)
	#print(sl_idx, "slatom_SASA", slatom_SASA, sl_coord, "cplx_SASA", cplx_sl_SASA, cplx_sl_coord)
	superligand_sasa_total = superligand_sasa_total + float(slatom_SASA)
	iidx = supligatom.GetIdx()
	#print(iidx)
	diff = float(slatom_SASA) - float(cplx_sl_SASA)
	#print("DIFF:", float(slatom_SASA) - float(cplx_sl_SASA))
	if (diff < 15.0) and  (float(cplx_sl_SASA) > 0.0): #flag
		mi  =  Chem.AtomPDBResidueInfo(orig_sl_info.GetName())
		mi.SetResidueName('TRP')
		mi.SetTempFactor(1.0)
		sl_atom.SetMonomerInfo(mi)
	else:
		mi  =  Chem.AtomPDBResidueInfo(orig_sl_info.GetName()) #If the name is not set, the PDB file will not be displayed correctly in PyMOL
		mi.SetResidueName('ASP')
		mi.SetTempFactor(0.0)
		sl_atom.SetMonomerInfo(mi)

	if cplx_sl_SASA != slatom_SASA:
		superligand_delta_sasa_total = superligand_delta_sasa_total + (float(slatom_SASA) - float(cplx_sl_SASA))

	elif slatom_SASA != 0: #atoms which have 0 SASA are buried inside the sl
		buried_sl_atom += 1 #but these atoms do not reside on the surface
		
print(superligand_sasa_total, 'superligand surface area')
print(superligand_delta_sasa_total, 'superligand buried area')

Chem.rdmolfiles.MolToPDBFile(superligand, "sl_" + id + "_exp_atoms_modified.pdb")

# Calculate final descriptors
fraction_sasa_change = superligand_delta_sasa_total/superligand_sasa_total #fsasa
# the fraction_sasa_change might not correlate with druggability, because it will be the same for a druggable and a nondruggable site of similar shapes
# but proportionally different sizes. Perhaps we could consider the total surface area of superligand buried inside the pocket instead? Done below!
not_buried_sasa = superligand_sasa_total - superligand_delta_sasa_total


csa = total_delta_sasa
print("csa", csa)
hsa_t = hydrophobic_sasa_total
hsa_t_r = hsa_t / csa
print("hsa_t", hsa_t)
hsa_plus_aromatic_t = aromatic_N_O_sasa + hsa_t
print("hsa_plus_aromatic_t", hsa_plus_aromatic_t)
hydrophilic_sasa_total = total_delta_sasa-hsa_t
print("hydrophilic_sasa_total", hydrophilic_sasa_total)
csa_vol_r = float(csa) / float(vol)
sa_vol_r = float(sl_surf) / float(vol)



psa_r = hydrophilic_sasa_total/total_delta_sasa
psa_don_acc_r = don_acc_sasa / float(total_delta_sasa)
hiaa = sum_hydrophobicity_index/float(num_res)
haa = num_apolar_res/num_res
asa = aromatic_sasa_total
chsa = charged_sasa_total 
asa_r = asa/csa
chsa_r = chsa/csa
sl_csa = superligand_delta_sasa_total
fr_buried_sl_atoms = buried_sl_atom / surface_atoms #old name: fr_exp_sl_atoms
fr_surf_sl_atoms = surface_atoms / superligand.GetNumAtoms()
no_bs_atoms = len(bs_atom_list)
sl_bs_r =  superligand.GetNumAtoms() / no_bs_atoms #

ph_t =  psa_don_acc_r + hsa_t_r # the sum of the polar and hydrophobic surface area ratios should be 1.

print('csa: ', csa, 'hsa_t: ', hsa_t, 'psa_r: ', psa_r, 'hiaa:', hiaa, 'haa:', haa, "chsa:", chsa, "f_C_atoms",str(number_c_atoms/no_bs_atoms), "C_atoms", number_c_atoms, "Num atoms", prot.GetNumAtoms())
print('fraction_sasa_change:', fraction_sasa_change, 'dsasa: ', not_buried_sasa)
print('asa: ', asa) #'asa_r: ', asa_r, 'chsa: ', chsa, 'chsa_r: ', chsa_r, 'psa_don_acc_r: ', psa_don_acc_r
print('surface_atoms', surface_atoms, 'buried_sl_atom', buried_sl_atom, 'total sl atoms', superligand.GetNumAtoms())

# Upload descriptors into new format table

#			id 	                csa 	            sl_csa 	               hsa 	                   exp_sl_sa 	                   fsasa 	                           asa 	                asa_r 	               chsa 	              chsa_r 	              psa_r 	                    aliphat_aromat_t 	                   fr_hpb_atoms 	                         c_ali_sa_r 	                       ali_sa_r 	                   binding_site_seq 	                         SpherocityIndex 	                                       Asphericity 	                                      Eccentricity 	              InertialShapeFactor 	                                                               RadiusOfGyration             	                        NPR2 	                                     NPR1 	                                      PMI3 	                                       PMI2 	                                    PM1 	              no_sl_atoms 	                   no_bs_atoms 	                fr_buried_sl_atoms 	               sl_bs_r 	               fr_surf_sl_atoms 	             vol 	             csa_vol_r 	               sa_vol_r 
#line =  str(id  + " \t "  + str(csa) + " \t " + str(sl_csa) + " \t " + str(hsa_t)  + " \t "  + str(not_buried_sasa) + " \t " + str(fraction_sasa_change) + " \t "  + str(chsa) + " \t "  + str(chsa_r)  + " \t "  + str(psa_don_acc_r)+ " \t "  + str(hsa_plus_aromatic_t) + " \t "  + str(number_c_atoms/bs_atom) + " \t "  + str(c_aliphatic_sasa/csa)  +  " \t " +str(all_aliphatic_sasa/csa)+ " \t " + str(bs_seq) + " \t " + str(Descriptors3D_values['SpherocityIndex']) + " \t " + str(Descriptors3D_values['Asphericity']) + " \t " + str(Descriptors3D_values['Eccentricity']) + " \t " + str(Descriptors3D_values['InertialShapeFactor']) + " \t " + str(Descriptors3D_values['RadiusOfGyration']) + " \t " + str(Descriptors3D_values['NPR2']) + " \t " + str(Descriptors3D_values['NPR1']) + " \t " + str(Descriptors3D_values['PMI3']) + " \t " + str(Descriptors3D_values['PMI2']) + " \t " + str(Descriptors3D_values['PMI1']) + " \t "+ str(no_sl_atoms) + " \t " + str(len(bs_atom_list)) + " \t " + str(fr_buried_sl_atoms) + " \t " + str(sl_bs_r) + " \t " + str(fr_surf_sl_atoms) +  " \t " + str(vol) + " \t " + str(csa_vol_r) + " \t " + str(sa_vol_r) + "\n")
#print(line)
#out.write(line)

#         ["id" ,   "csa" ,       "sl_csa" ,   "hsa" ,     "exp_sl_sa" ,           "fsasa" ,  "chsa" ,   "chsa_r" ,   "psa_r" ,   "aliphat_aromat_t" ,   "fr_hpb_atoms" ,   "c_ali_sa_r" ,   "ali_sa_r" ,   "binding_site_seq" ,   "SpherocityIndex" ,   "Asphericity" ,   "Eccentricity" ,   "InertialShapeFactor" ,   "RadiusOfGyration" ,   "NPR2" ,   "NPR1" ,   "PMI3" ,   "PMI2" ,   "PM1" ,   "no_sl_atoms" ,   "no_bs_atoms" ,   "fr_buried_sl_atoms" ,   "sl_bs_r" ,   "fr_surf_sl_atoms" ,   "vol" ,   "csa_vol_r" ,   "sa_vol_r" ] )
line = [str(id) ,str(csa) ,  str(sl_csa) ,  str(hsa_t) ,str(not_buried_sasa) ,  str(fraction_sasa_change) ,str(chsa) ,str(chsa_r) ,str(psa_don_acc_r),str(hsa_plus_aromatic_t) ,str(number_c_atoms/bs_atom) ,str(c_aliphatic_sasa/csa)  , str(all_aliphatic_sasa/csa),  str(bs_seq) ,  str(Descriptors3D_values['SpherocityIndex']) ,  str(Descriptors3D_values['Asphericity']) ,  str(Descriptors3D_values['Eccentricity']) ,  str(Descriptors3D_values['InertialShapeFactor']) ,  str(Descriptors3D_values['RadiusOfGyration']) ,  str(Descriptors3D_values['NPR2']) ,  str(Descriptors3D_values['NPR1']) ,  str(Descriptors3D_values['PMI3']) ,  str(Descriptors3D_values['PMI2']) ,  str(Descriptors3D_values['PMI1']) , str(no_sl_atoms) ,  str(len(bs_atom_list)) ,  str(fr_buried_sl_atoms) ,  str(sl_bs_r) ,  str(fr_surf_sl_atoms) ,  str(vol) ,  str(csa_vol_r) ,  str(sa_vol_r) ]


with open("../../../descriptor_values.csv", "a") as out:

	writer = csv.writer(out, delimiter = ",")
	print(line)
	#out.write("id \t csa \t sl_csa \t hsa \t exp_sl_sa \t fsasa \t asa \t asa_r \t chsa \t chsa_r \t psa_r \t aliphat_aromat_t \t fr_hpb_atoms \t c_ali_sa_r \t ali_sa_r \t  binding_site_seq \t  SpherocityIndex \t  Asphericity \t  Eccentricity \t  InertialShapeFactor \t  RadiusOfGyration \t  NPR2 \t  NPR1 \t  PMI3 \t  PMI2 \t  PM1 \t  no_sl_atoms \t  no_bs_atoms \t  fr_buried_sl_atoms \t  sl_bs_r \t  fr_surf_sl_atoms \t  vol \t  csa_vol_r \t  sa_vol_r \n")
	writer.writerow(line) 

out.close()
#command = "INSERT INTO " + dt + " (id,csa,sl_csa,hsa,psa_r_old,hiaa,haa,exp_sl_sa,fsasa,asa,asa_r,chsa,chsa_r,psa_r,aliphat_aromat_t,fr_hpb_atoms,c_ali_sa_r,ali_sa_r, binding_site_seq, SpherocityIndex, Asphericity, Eccentricity, InertialShapeFactor, RadiusOfGyration, NPR2, NPR1, PMI3, PMI2, PM1, dsasa_2, dsasa_3, no_sl_atoms, no_bs_atoms, fr_buried_sl_atoms, sl_bs_r, fr_surf_sl_atoms, proberad, vol, csa_vol_r, sa_vol_r, don_acc_sasa, ph_t) values ( '" + id  + "',"  + str(csa) + "," + str(sl_csa) + "," + str(hsa_t) + "," + str(psa_r)  + "," + str(hiaa) + "," + str(haa) + "," + str(not_buried_sasa) + "," + str(fraction_sasa_change) + ","  + str(asa) + ","  + str(asa_r) + ","  + str(chsa) + ","  + str(chsa_r)  + ","  + str(psa_don_acc_r)+ ","  + str(hsa_plus_aromatic_t) + ","  + str(number_c_atoms/bs_atom) + ","  + str(c_aliphatic_sasa/csa)  +  "," +str(all_aliphatic_sasa/csa)+ "," + "'%s'" % bs_seq + "," + str(Descriptors3D_values['SpherocityIndex']) + "," + str(Descriptors3D_values['Asphericity']) + "," + str(Descriptors3D_values['Eccentricity']) + "," + str(Descriptors3D_values['InertialShapeFactor']) + "," + str(Descriptors3D_values['RadiusOfGyration']) + "," + str(Descriptors3D_values['NPR2']) + "," + str(Descriptors3D_values['NPR1']) + "," + str(Descriptors3D_values['PMI3']) + "," + str(Descriptors3D_values['PMI2']) + "," + str(Descriptors3D_values['PMI1']) + "," +  str(dsasa_2) + "," + str(dsasa_3) + "," + str(no_sl_atoms) + "," + str(len(bs_atom_list)) + "," + str(fr_buried_sl_atoms) + "," + str(sl_bs_r) + "," + str(fr_surf_sl_atoms) + "," + str(opts.probeRadius) + "," + str(vol) + "," + str(csa_vol_r) + "," + str(sa_vol_r) + "," + str(don_acc_sasa) + "," + str(ph_t) + ")"


#print(command)
#cursor.execute(command)
#cursor.close ()
#conn.commit()
#conn.close ()
