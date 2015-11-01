#!/usr/bin/env python
import sys
import os
import re
import subprocess
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

def sortIt_length(file, verbose_level=0):
	"""Sorts clusters by seq length"""
	outFile = re.sub(r"(.*)\..*", r"\1_Sl.fa", file) # Make the outfile name by cutting off the extension of the infile name, and adding "_S1.fa"
	usearch_cline = "%s -sortbylength %s -fastaout %s" %(Usearch, file, outFile)
	#usearch_cline = "%s -sortbylength %s -output %s" %(Usearch, file, outFile) # Usearch 7 
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		#print '\n**Usearch-sorting output on', file, '**\n'
		#print err
		log.write('\n**Usearch-sorting output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile # having it spit out the outfile name, if necessary, so that it can be used to call downstream stuff and avoid complicated glob.globbing

def clusterIt(file, clustID, round, verbose_level=0):
	"""The clustering step, using the clustID value"""
	outFile = re.sub(r"(.*)\.fa", r"\1C%s_%s.fa" %(round, clustID), file) # The rs indicate "raw" and thus python's escaping gets turned off
	outClustFile = re.sub(r"(.*)\.fa", r"\1clusts%s\.uc" %(round), file)	
	if round == 1:
		usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -sortedby other -centroids %s -uc %s -sizeout" % (Usearch, file, clustID, outFile, outClustFile) 
		#usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -usersort -consout %s -uc %s -sizeout" % (Usearch, file, clustID, outFile, outClustFile) # Usearch 7   
        # Can add in "-cons_truncate" to the usearch call, if the primer removal isn't effective, but note some problems with partial sequences results.
	elif round > 1:
		usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -sortedby other -centroids %s -uc %s -sizein -sizeout" % (Usearch, file, clustID, outFile, outClustFile)
		#usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -usersort -consout %s -uc %s -sizein -sizeout" % (Usearch, file, clustID, outFile, outClustFile) # Usearch 7	
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		#print '\n**Usearch-clustering output on', file, '**\n'
		#print err
		log.write('\n**Usearch-clustering output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile
	
# This function looks for PCR chimeras -- those formed within a single amplicon pool
def deChimeIt(file, round, verbose_level=0):
	"""Chimera remover. The abskew parameter is hardcoded currently (UCHIME default for it is 2.0)"""
	outFile = re.sub(r"(.*)\.fa", r"\1dCh%s.fa" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
	usearch_cline = "%s -uchime_denovo %s -abskew 1.9 -nonchimeras %s" % (Usearch, file, outFile)
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level in [1, 2]:
		#print '\n**Uchime output on', file, '**\n'
		#print err
		log.write('\n**Uchime output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile
	
def sortIt_size(file, thresh, round, verbose_level=0):
	"""Sorts clusters by size, and removes those that are smaller than a particular size
	(sent as thresh -- ie, sizeThreshold or sizeThreshold2). 
    "round" is used to annotate the outfile name with S1, S2, etc. depending on which sort this is"""
	outFile = re.sub(r"(.*)\.fa", r"\1Ss%s.fa" %(round), file)
	usearch_cline = "%s -sortbysize %s -fastaout %s -minsize %f" %(Usearch, file, outFile, thresh)
	#usearch_cline = "%s -sortbysize %s -output %s -minsize %f" %(Usearch, file, outFile, thresh) # Usearch 7
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		#print '\n**Usearch-sorting output on', file, '**\n'
		#print err
		log.write('\n**Usearch-sorting output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile

def ClusterDechimera(annotd_seqs_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2, Multiplex_per_barcode = False, verbose_level = 0):
	sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
	locusCounts = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus", Multiplex_perBC_flag = Multiplex_per_barcode) 

	## Split the locus files by barcode, and cluster each of the resulting single locus/barcode files
	sys.stderr.write('Clustering/dechimera-izing seqs...\n')
	log.write('#Sequence clustering/dechimera-izing#\n')
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
			if Multiplex_per_barcode:	
				bcodeCounts = SplitBy(annotd_seqs_file = each_folder + ".fa", split_by = "barcode", Multiplex_perBC_flag = Multiplex_per_barcode)
				all_folders_bcodes = bcodeCounts.keys()
			else:
				taxonCounts = SplitBy(annotd_seqs_file = each_folder + ".fa", split_by = "taxon", Multiplex_perBC_flag = Multiplex_per_barcode)					
				all_folders_bcodes = taxonCounts.keys()

			for bcode_folder in all_folders_bcodes:
				if verbose_level in [1,2]:
					log.write("Working on " + bcode_folder + '\n')
				os.chdir(bcode_folder)
				
				if verbose_level in [1,2]:
					log.write("\tFirst clustering\n")
					log.write("\nAttempting to sort: " + bcode_folder + ".fa\n")

				sorted_length = sortIt_length(file = bcode_folder + ".fa", verbose_level = verbose_level)
				clustered1 = clusterIt(file = sorted_length, clustID = clustID, round = 1, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tFirst chimera slaying expedition\n")
				deChimered1 = deChimeIt(file = clustered1, round = 1, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tSecond clustering\n")
				sorted_size1 = sortIt_size(file = deChimered1, thresh = sizeThreshold, round = 1, verbose_level = verbose_level)
				clustered2 = clusterIt(file = sorted_size1, clustID = clustID2, round = 2, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tSecond chimera slaying expedition\n")
				deChimered2 = deChimeIt(file = clustered2, round = 2, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tThird clustering\n")
				sorted_size2 = sortIt_size(file = deChimered2, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
				clustered3 = clusterIt(file = sorted_size2, clustID = clustID3, round = 3, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tThird chimera slaying expedition\n")
				deChimered3 = deChimeIt(file = clustered3, round = 3, verbose_level = verbose_level)

				# Why are we sorting again? I guess this gives us the chance to remove clusters smaller than sizeThreshold2
				sorted_size3 = sortIt_size(file = deChimered3, thresh = sizeThreshold2, round = 3, verbose_level = verbose_level)

				try:
					clustered_seq_file = parse_fasta(sorted_size3)
					for each_seq in clustered_seq_file:
						taxon_name = str(each_seq.id).split('|')[0].split('=')[-1] # for example, get C_dia_5316 from centroid=centroid=C_dia_5316|ApP|C|BC02|_p0/158510/ccs;ee=1.9;;seqs=6;seqs=18;size=27;
						
						try:
							LocusTaxonCountDict_clustd[taxon_name, each_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
							# The locus names are the same as each_folder
						except:
							LocusTaxonCountDict_clustd[taxon_name, each_folder] = 1		
				except:
					if verbose_level in [1,2]:
						log.write(str(sorted_size3) + 'is an empty file\n')

				os.chdir("..") # To get out of the current barcode folder and ready for the next one
		os.chdir("..") # To get out of the current locus folder and ready for the next one
	log.write('\t...done\n\n')	

	return
ppp_location = os.path.dirname(os.path.abspath( __file__ ))
Usearch = ppp_location + '/' + 'Dependencies/usearch8.1.1756'

log = open('log.txt', 'w')
annotated_file = sys.argv[1]
#split_type = sys.argv[2]

ClusterDechimera(annotated_file, 0.997, 0.995, 0.990, 1, 4)

"""
try:
	Multiplex_perBC_flag = sys.argv[3]
	if Multiplex_perBC_flag == '-M':
		SplitBy(annotated_file, split_type, Multiplex_perBC_flag=True)
	else:
		print 'Error'
except:
	SplitBy(annotated_file, split_type, Multiplex_perBC_flag=False)
"""












