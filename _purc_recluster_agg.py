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

import sys
import os
import re
import subprocess
import glob
import shutil
import time
import datetime
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline


usage = """

Use this script to recluster the alleles/homeologs from a previous PURC run. 

Usage: ./purc_recluster_agg.py annotated_file output_folder clustID1 sizeThreshold
Example: ./purc_recluster_agg.py purc_run_3_annotated.fa Run2 0.99 4

Note: 
(1) clustID: The similarity criterion for the first, second and third USEARCH clustering
(2) sizeThreshold: The min. number of sequences/cluster necessary for that cluster to be retained (set to 2 to remove singletons, 3 to remove singletons and doubles, etc)

	"""

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
	LocusTaxonCountDict_unclustd = {}
	splits_list = []

	for each_seq in unsplit_seq:
		#finding the identifing annotation for the split of interest.
		# e.g., BC01, BC02, ... or GAP, PGI, ...

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
	
		if split_by == "locus":
			try:
				LocusTaxonCountDict_unclustd[str(each_seq.id).split('|')[0], str(each_seq.id).split('|')[1]] += 1  # {('C_dia_5316', 'ApP'): 28} for example
				# The locus names are the same as each_folder
			except:
				LocusTaxonCountDict_unclustd[str(each_seq.id).split('|')[0], str(each_seq.id).split('|')[1]] = 1		
		
		os.chdir('..')
	if split_by == "locus":
		return splits_count_dic, LocusTaxonCountDict_unclustd #as {'BC01': 150, 'BC02': 156} for example
	else:
		return splits_count_dic

def align_and_consensus(inputfile, output_prefix):
	output_alignment = output_prefix + '_aligned.fa'
	output_consensus = output_prefix + '_consensus.fa'
	muscle_cline = '%s -in %s -out %s' % (Muscle, inputfile, output_alignment)
	process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	#muscle_cline = MuscleCommandline(input=inputfile, out=output_alignment) 
	#stdout, stderr = muscle_cline()

	alignment = AlignIO.read(output_alignment, 'fasta')
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.gap_consensus(ambiguous='N')
	output = open(output_consensus, 'w')
	output.write('>' + output_prefix + '_consensus' + '\n' + str(consensus).replace('-','') + '\n')
	return

def clusterIt_agg(file, clustID, sizeThreshold, verbose_level=2):
	#clustering
	round = '1'
	outClustFile = re.sub(r"(.*)\.fa", r"\1C%s_%s.txt" %(round, clustID), file) # The rs indicate "raw" and thus python's escaping gets turned off	
	outClustFile_seq = re.sub(r"(.*)\.fa", r"\1C%s_%s.fa" %(round, clustID), file) # The rs indicate "raw" and thus python's escaping gets turned off		
	usearch_cline = "%s -cluster_agg %s -id %f -clusterout %s -linkage max -treeout tree.tre" % (Usearch, file, clustID, outClustFile) 

	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		#print '\n**Usearch-clustering output on', file, '**\n'
		#print err
		log.write('\n**Usearch-clustering output on ' + str(file) + '**\n')
		log.write(str(err))
	
	#parse the outClustFile to get sequences
	outClust = open(outClustFile, 'rU')
	ClustToSeq_dict = {}
	for line in outClust:
		line = line.strip('\n')
		cluster = line.split('\t')[0]
		sequence = line.split('\t')[1]
		try:
			ClustToSeq_dict[cluster].append(sequence)
		except:
			ClustToSeq_dict[cluster] = [sequence]
	#print ClustToSeq_dict
	SeqDict = SeqIO.index(file, 'fasta')
	for cluster in ClustToSeq_dict:
		if len(ClustToSeq_dict[cluster]) > int(sizeThreshold):
			cluster_seq_file = open(cluster, 'w')
			for seq in ClustToSeq_dict[cluster]:
				cluster_seq_file.write('>' + seq + '\n' + str(SeqDict[seq].seq) + '\n')
			cluster_seq_file.close()
			align_and_consensus(cluster, file.split('.fa')[0])

def ClusterDechimera(annotd_seqs_file, clustID, sizeThreshold, verbose_level = 1): # M_p_barcode was set to FALSE, whcih led to splitting by taxon instead of barcode
	sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
	locusCounts, LocusTaxonCountDict_unclustd = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus") 

	## Split the locus files by barcode, and cluster each of the resulting single locus/barcode files
	all_folders_loci = locusCounts.keys() # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
		# folders) and they correspond to the counts for each.
	LocusTaxonCountDict_clustd = {}

	for each_folder in all_folders_loci: 
		os.chdir(each_folder)
		sys.stderr.write('\nWorking on: ' + each_folder + '...\n')
		if verbose_level in [1,2]:
			log.write('\nWorking on ' + str(each_folder) + ' ...\n')
		# Open the file with annotated sequences for that locus.This is a little awkward, but the file name is the same as the folder name (which is "each_folder" currently)
		# with the addition of ".fa". This is set as the file handle in SplitsBy, and if changed there, needs to be changed here too
		AnnodDict = SeqIO.index(each_folder + ".fa", 'fasta') 
		if len(AnnodDict) > 0: # ie, the file is not empty
			taxonCounts = SplitBy(annotd_seqs_file = each_folder + ".fa", split_by = "taxon")					
			all_folders_bcodes = taxonCounts.keys()

			for bcode_folder in all_folders_bcodes:
				if verbose_level in [1,2]:
					log.write("Working on " + bcode_folder + '\n')
				os.chdir(bcode_folder)
				clusterIt_agg(bcode_folder + ".fa", 0.99, 4, verbose_level=2)


				os.chdir("..") # To get out of the current barcode folder and ready for the next one
		os.chdir("..") # To get out of the current locus folder and ready for the next one
	log.write('\t...done\n\n')	

	## Put all the sequences together ##
	sys.stderr.write('\rPutting all the sequences together......\n\n')
	for each_folder in all_folders_loci: # Looping through each of the locus folders
		outputfile_name = str(each_folder) + '_clustered.fa' # "each_folder" is also the name of the locus
		outputfile = open(outputfile_name, 'w')
		os.chdir(each_folder)
		bcodesForThisLocus = glob.glob("*")
		#print "glob thinks there are these barcode folders present: ", bcodesForThisLocus, "\n\n"
		for bcode_folder in bcodesForThisLocus: # have to go into each barcode folder in each locus folder
			if os.path.isdir(bcode_folder): # the glob might have found some files as well as folders	
				os.chdir(bcode_folder)
				files_to_add = glob.glob("*consensus.fa")
				for file in files_to_add: 
					shutil.copyfileobj(open(file,'rb'), outputfile) #Add each file to the final output		
				os.chdir('..')
		os.chdir('..')
		outputfile.close()

	return LocusTaxonCountDict_clustd, LocusTaxonCountDict_unclustd

def muscleIt(file, verbose_level=0):
	"""Aligns the sequences using MUSCLE"""
	outFile = re.sub(r"(.*)\..*", r"\1.afa", file) # The rs indicate "raw" and thus python's escaping gets turned off
	muscle_cline = '%s -in %s -out %s' % (Muscle, file, outFile)
	process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		#print '\n**Muscle output on', file, '**\n'
		#print err
		log.write('\n**Muscle output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile

###### RUN ######	
if len(sys.argv) < 2:
	sys.exit(usage)

annotated_file = sys.argv[1]
masterFolder = sys.argv[2]
clustID = float(sys.argv[3])
sizeThreshold = int(sys.argv[4])

purc_location = os.path.dirname(os.path.abspath( __file__ ))
Usearch = purc_location + '/' + 'Dependencies/usearch8.1.1756'
Muscle = purc_location + '/' + 'Dependencies/muscle3.8.31'

## Make output folder ##
if os.path.exists(masterFolder): # overwrite existing folder
	shutil.rmtree(masterFolder)
os.makedirs(masterFolder)
os.chdir(masterFolder)

ts = time.time()
time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
log_file = 'purc_log_' + time_stamp + '.txt'
log = open(log_file, 'w')
log.write(logo + '\n')
log.write('purc_recluster_agg.py ' + annotated_file + ' ' + masterFolder + ' ' + str(clustID) + ' ' + str(clustID2) + ' ' + str(clustID3) + ' ' + str(sizeThreshold) + ' ' + str(sizeThreshold2) + '\n\n')

## Recluster and redechimera ##
ClusterDechimera('../'+annotated_file, clustID, sizeThreshold)

'''
taxon_list = []
locus_list = []
for taxon_locus in LocusTaxonCountDict_clustd.keys():
	taxon_list.append(taxon_locus[0])
	locus_list.append(taxon_locus[1])
taxon_list = set(taxon_list)
locus_list = set(locus_list)

## Count ##
count_output = open('purc_cluster_counts.xls', 'w')

log.write("\n**Raw reads per accession per locus**\n")
count_output.write('\n**Raw reads per accession per locus**\n')
count_output.write('\t' + '\t'.join(locus_list) + '\n')
log.write('\t' + '\t'.join(locus_list) + '\n')
for each_taxon in set(taxon_list): 
	count_output.write(each_taxon + '\t')
	log.write(each_taxon + '\t')
	for each_locus in locus_list:
		try:
			count_output.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
			log.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
		except:
			count_output.write('0\t')
			log.write('0\t')
	count_output.write('\n')
	log.write('\n')

count_output.write('\n**Final clustered sequences per accession per locus**\n')
log.write('\n**Final clustered sequences per accession per locus**\n')
count_output.write('\t' + '\t'.join(locus_list) + '\n')	
log.write('\t' + '\t'.join(locus_list) + '\n')	
for each_taxon in set(taxon_list):
	count_output.write(each_taxon + '\t')		
	log.write(each_taxon + '\t')		
	for each_locus in locus_list:
		try:
			count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
			log.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
		except:
			count_output.write('0\t') 
			log.write('0\t') 
	count_output.write('\n')		
	log.write('\n')		

count_output.write('\n**Allele/copy/cluster/whatever count by locus**\n')	
log.write('\n**Allele/copy/cluster/whatever count by locus**\n')	

for each_locus in locus_list:
	file_name = str(each_locus) + '_clustered.txt'
	try: 
		seq_no = len(parse_fasta(file_name))
		count_output.write(str(each_locus) + '\t' + str(seq_no) + '\n')
		log.write(str(each_locus) + '\t' + str(seq_no) + '\n')
	except:
		count_output.write(str(each_locus) + '\t0\n')			
		log.write(str(each_locus) + '\t0\n')			

## Clean-up the sequence names ##
sys.stderr.write("Cleaning up the file names\n")
fastas = glob.glob("*_clustered.txt")
for file in fastas:
	#log.write("Cleaning up the sequence names in " + str(file) + "\n")
	fasta_cleaned = open(str(file).replace(".txt", "_renamed.fa"), 'w') # Make a new file with .fa instead of .txt
	parsed = parse_fasta(file)
	for seq in parsed:
		seq.id = re.sub(r"seqs=\d*", r"", seq.id)
		seq.id = re.sub(r"ccs.ee=[\d\.]*", r"", seq.id)
		seq.id = seq.id.replace("centroid=", "").replace(";", "").replace("/","_")
		fasta_cleaned.write(str('>' + seq.id + "\n" + seq.seq + "\n"))
	fasta_cleaned.close()

## Aligning the sequences ##
fastas = glob.glob("*_renamed.fa")
for file in fastas:
	sys.stderr.write("Aligning " + file + "\n")
	log.write("Aligning " + file + "\n")
	outFile = muscleIt(file)
	





'''




