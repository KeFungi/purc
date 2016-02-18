#!/usr/bin/env python

'''This script takes each locus alignment output for each PURC clustering regimes and uses MUSCLE to profile-profile
align them together. Output is one alignment per locus, which includes all sequences for that locus from each of 
the regimes'''

import glob
from Bio import SeqIO
import os

def parse_fasta(infile):
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq]

#workDir is the directory containing the subdirectories foe each of the analyses being compared
workDir = "/Users/carlrothfels/Desktop/PURC_analysesForPaper/6_runsWchimeraCountsOutput"
saveToDir = "/Users/carlrothfels/Desktop/PURC_analysesForPaper/6_runsWchimeraCountsOutput/replicabilityTests"

runs = ["R2"] # Presumably best to do each run (each pacbio dataset) independently to avoid huge alignments, 
# but could do multiple at once if desired.

regimes = ["a", "c", "e"] # However your clustering regimes are labelled
regimesToAdd = []
for regime in regimes: # Adding the "Str_a" etc. stringent regimes (this is only necessary because of the way I labelled my regimes)
	regimesToAdd.append("Str_" + regime)
regimes = regimes + regimesToAdd # Adding the stringent regimes to the list of regimes

for run in runs:
	
	for regime in regimes:
		directory = workDir + "/" + run + regime + "/" # E.g., "/Users/carlrothfels/Desktop/PURC_analysesForPaper/6_runsWchimeraCountsOutput/R2D/"
		os.chdir(directory)
		
		file_list = glob.glob('*reconsensus.afa') #In the output from PURC, there are two .afa files for each locus -- the raw clusters, and the new "reclustered" one

		for file in file_list:
			print "Looking in file: ", file
			os.chdir(directory) # getting back into the directory w the data
			ID_to_add = run + regime # Will annotate the sequences with this information, so that it's clear which run each seq. came from
			locus = file.split('_')[0] # Splitting file name at understores and taking first element. Works on files named like "PGI_reconsensus.afa"
			
			sequences = parse_fasta(file)
			os.chdir(saveToDir) # getting back into the main directory to store the output files
			
			# Need the script to work with names like:
			# >G_oya_6399_Cluster0;size=39;   (most of them) and:
			# >xCyst_7974_Cluster2;size=27;   (only this one doesn't have both the genus initial and the abbreviated specific epithet)

			for seq in sequences:
				#accession = str(seq.id).split('|')[0].split('_')[-1]
				if seq.id[0] == "x": # xCysto is the only seq name to start with an "x"
					accession = str(seq.id).split('|')[0].split('_')[1] #xCystocarpium seqs are labelled ">xCysto_7974" so only have one underscore
				else:
					accession = str(seq.id).split('|')[0].split('_')[2] #Trying to split by pipes and underscores and then take the third element					
					
				file_to_save = open(locus + '_repbiltyTest.fa', 'a')			
				file_to_save.write('>' + ID_to_add + '|' + str(seq.id) + '\n' + str(seq.seq) + '\n')
				file_to_save.close()
				## need to check here if there is already an alignment for this locus, and if there is, do a profile-profile alignment with it


