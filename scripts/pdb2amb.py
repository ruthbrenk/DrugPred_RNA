#!/usr/bin/python

# format pdb file with hydrogens to amber format for DOCK
# amber file must not contain sequence identifier => will be deleted by this script
import string,os, sys

import sys, string

in_file = open(sys.argv[1], 'r')
out_file = open(sys.argv[2], 'w')


for i in in_file.readlines():
	if i[:4] =="ATOM":


		#print 'test'
		new_line = i[:-1]
		#ALL
		new_line = new_line.replace('H   ', 'HN  ')
		#ARG
		new_line = new_line.replace(' HE ', ' HNE')
		new_line = new_line.replace('1HH1 ', ' HN11')
		new_line = new_line.replace('1HH2 ', ' HN12')
		new_line = new_line.replace('2HH1 ', ' HN21')
		new_line = new_line.replace('2HH2 ', ' HN22')
		new_line = new_line.replace('HH11 ', ' HN11')
		new_line = new_line.replace('HH12 ', ' HN12')
		new_line = new_line.replace('HH21 ', ' HN21')
		new_line = new_line.replace('HH22 ', ' HN22')
		#ASN
		new_line = new_line.replace('1HD2 ', ' HND1')
		new_line = new_line.replace('2HD2 ', ' HND2')
		new_line = new_line.replace('HD22 ', ' HND1')
		new_line = new_line.replace('HD21 ', ' HND2')
		#CYS
		new_line = new_line.replace('HG  CYS', 'HSG CYS')
		if new_line.find('SER') != -1:
			print (new_line)


		#SER
		new_line = new_line.replace('HG  SER', 'HOG SER')
		#GLN
		new_line = new_line.replace('1HE2 ', ' HNE1')
		new_line = new_line.replace('2HE2 ', ' HNE2')
		new_line = new_line.replace('HE22 ', ' HNE1')
		new_line = new_line.replace('HE21 ', ' HNE2')
		#HIS
		new_line = new_line.replace('HD1', 'HND')
		#LYS
		new_line = new_line.replace('HZ1 ', 'HNZ1')
		new_line = new_line.replace('HZ2 ', 'HNZ2')
		new_line = new_line.replace('HZ3 ', 'HNZ3')
		#TRP
		new_line = new_line.replace('HE1', 'HNE')
		#TYR
		new_line = new_line.replace('HH ', 'HOH')
		#THR
		new_line = new_line.replace('HG1', 'HOG')

		#Bases
		new_line = new_line.replace('H2 1', ' H21')
		new_line = new_line.replace('H2 2', ' H22')
		new_line = new_line.replace(' HN1', ' H1 ')
		new_line = new_line.replace('H6 1', ' H61')
		new_line = new_line.replace('H6 2', ' H62')
		new_line = new_line.replace(' HN3', ' H3 ')
		new_line = new_line.replace('H4 1', ' H41')
		new_line = new_line.replace('H4 2', ' H42')
		new_line = new_line.replace('H5\'1', ' H5\'')
		new_line = new_line.replace('H5\'2', ' H5\'')
		new_line = new_line.replace('HO2\'', ' HO2')

		#out_file.write(new_line + '\n')
		#has to be here because otherwise we will lose information for ARG residues



		new_line = new_line[:12] + new_line[13:17] + ' ' + new_line[17:21] + ' ' + new_line[22:55] + '\n'



		out_file.write(new_line)
	else:
		out_file.write(i)
	
