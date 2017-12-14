#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:30:55 2017. Updated Wed Dec 13, 2017

@author: jonathangreenhalgh
"""

#Script to model sequence of a protein
#Required Inputs:
#   (root_name).ali:an alignment or sequence in PIR format. 
#   root_name:      providing a root_name as a string will help the program 
#                   find the alignment file and name the outputs to match. 
#                   This should match the name of the alignment file (
#                   without an extension)
#   'pdb_95.pir':   A file that contains non-redundant entries in the pdb
#   '
#Outputs:
#   pdb_(root_name).bin:  binary database output
#   build_(root_name)_profile.prf
#   build_(root_name)_profile.ali:  Profile alignment file in PIR format. Used
#                                   to evaluate identity, coverage etc of 
#                                   sequences in pdb_95.pir
#
#   build_(root_name)_profile_parameters.csv:
#       csv file containing the % sequence identity, coverage, and E scores of 
#       potential template sequences in the PDB.



from modeller import *
    
log.verbose()
env = environ()

#Root_filename
#if alignment file is named root.ali, root_name=root

root_name='model_protein' #This is a tag that will be incorporated in the filenames. 

#-- Prepare the input files

#-- Read in the sequence database
sdb = sequence_db(env)
sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 8000), clean_sequences=True)

#-- Write the sequence database in binary form
sdb.write(seq_database_file='pdb_'+root_name+'.bin', seq_database_format='BINARY',
          chains_list='ALL')

#-- Now, read in the binary database
sdb.read(seq_database_file='pdb_'+root_name+'.bin', seq_database_format='BINARY',
         chains_list='ALL')

#-- Read in the target sequence/alignment
aln = alignment(env)
aln.append(file=root_name+'.ali', alignment_format='PIR', align_codes='ALL') 

#-- Convert the input sequence/alignment into
#   profile format
prf = aln.to_profile()

#-- Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

#-- Write out the profile in text format
prf.write(file='build_'+root_name+'_profile.prf', profile_format='TEXT')

#-- Convert the profile back to alignment format
aln = prf.to_alignment()

#-- Write out the alignment file
aln.write(file='build_'+root_name+'_profile.ali', alignment_format='PIR')

#-- Read the prf file, and keep all numbers except the sequence of the aligned portion

short_prf=open('build_'+root_name+'_profile.prf','r')

short_prf_list=[]                   #make an empty list to populate
for line in short_prf:
    line=line.split(' ')            #split the line by spaces
    line=list(filter(None,line))    #remove all the empty entries in the list
    line=line[:-1]                  #remove the sequence of the aligned portion
    short_prf_list.append(line)     #append line to the list
    
short_prf_list=short_prf_list[6:]   #Cut out introductory information

for line in short_prf_list:
    line[0]=int(line[0])            #change the Number to an integer
    for j in range(3,9):            #change the other lengths to an integer
        line[j]=int(line[j])   
    line[9]=float(line[9])          #change the % id to a float
    line[10]=float(line[10])        #change the E value to a float

#Save the list in a more readable format
info_table=[['Number','PDB code','S/X','', '','','','','','Alignment length','% ID','E']]

for row in short_prf_list:
    info_table.append(row)
    
import numpy
numpy.savetxt('build_'+root_name+'_profile_parameters.csv',info_table,fmt='%s',delimiter=',')
