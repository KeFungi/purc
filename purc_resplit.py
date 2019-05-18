#!/usr/bin/env python

logo = """
-------------------------------------------------------------
|                            PURC                           |
|        Pipeline for Untangling Reticulate Complexes       |
|                        version 1.0                        | 
|            https://bitbucket.org/crothfels/ppp            |
|															|
|                 Fay-Wei Li & Carl J Rothfels              |
|           see purc.py script for more information         |
-------------------------------------------------------------
""" 

usage = """

Use this script to split an annotated fasta file based on taxon, barcode, locus, or group. 

Usage: ./purc_resplit.py annoated_file split_type
Example: ./purc_resplit.py purc_run_3_annotated.fa barcode

Note: 
(1) Possible split_type: barcode, locus, taxon, group
(2) If one barcode contains multiple samples, then pass the -M argument. e.g. ./purc_resplit.py purc_run_3_annotated.fa barcode -M  

"""
	
import sys
import os
from Bio import SeqIO

def parse_fasta(infile):
	"""Reads in a fasta, returns a list of biopython seq objects"""
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq] 

def SplitBy(annotd_seqs_file, split_by = "locus-taxon", Multiplex_perBC_flag=False): 
	"""Uses the annotated sequences to split sequences into different files based on splits_list 
	(could be by barcode or by locus, etc); returns a dictionary of seq counts for each subgroup"""
	
	#annotd_seqs = open(annotd_seqs_file, 'rU')
	unsplit_seq = parse_fasta(annotd_seqs_file)	
	splits_file_dic = {}
	splits_count_dic = {}
	splits_list = []

	for each_seq in unsplit_seq:
		#finding the identifing annotation for the split of interest.
		# e.g., BC01, BC02, ... or gapCp, PGIC, ...

		if split_by == "taxon":
			split = str(each_seq.id).split('|')[0]
		elif split_by == "locus":			
			split = str(each_seq.id).split('|')[1]
		elif split_by == "taxon-locus":
			split = str(each_seq.id).split('|')[0] + "_" + str(each_seq.id).split('|')[1]
		elif split_by == "locus-taxon": # same as above, but folders/files labeled with locus name before taxon name
			split = str(each_seq.id).split('|')[1] + "_" + str(each_seq.id).split('|')[0]
		elif split_by == "barcode" and Multiplex_perBC_flag:
			split = str(each_seq.id).split('|')[3]
		elif split_by == "group" and Multiplex_perBC_flag: # where group is the let set to identify the taxonomic groups, e.g, A, B, C, ...
			split = str(each_seq.id).split('|')[2]			
		elif split_by == "barcode" and not Multiplex_perBC_flag:
			split = str(each_seq.id).split('|')[2]
		else:
			sys.stderr.write('Error: Attempting to split sequences by an inappropriate (absent?) category.\n')

		if split not in splits_list: # add split identifier to the list of splits
			splits_list.append(split)
			os.makedirs(split)
			os.chdir(split)		
			seq_file = split + '.fa'
			file_handle = open(seq_file, 'w') #FWL only have to do this once, for the first time that split occurs?
			splits_file_dic[split] = file_handle # e.g., {BC01:BC01.fa, BC02:BC02.fa, ...}
			# os.chdir('..')
		else:
			os.chdir(split)	
	
		file_handle = splits_file_dic[split] #use that split as the key to find the corresponding output file
		file_handle.write('>' + str(each_seq.id) + '\n' + str(each_seq.seq) + '\n') #write the seq
		
		try:
			splits_count_dic[split] += 1
		except:
			splits_count_dic[split] = 1		
		
		os.chdir('..')
	# FWL - do we not have to close each of the files we've been writing to?
	return splits_count_dic #as {'BC01': 150, 'BC02': 156} for example



if len(sys.argv) < 3:
	sys.exit(usage)

annotated_file = sys.argv[1]
split_type = sys.argv[2]

try:
	Multiplex_perBC_flag = sys.argv[3]
	if Multiplex_perBC_flag == '-M':
		SplitBy(annotated_file, split_type, Multiplex_perBC_flag=True)
	else:
		print 'Error: unknown argument'
except:
	SplitBy(annotated_file, split_type, Multiplex_perBC_flag=False)













