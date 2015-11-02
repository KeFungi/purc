#!/usr/bin/env python
import sys
import os
import re
import subprocess
import glob
import shutil
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

	## Put all the sequences together ##
	sys.stderr.write('\rPutting all the sequences together......\n\n')
	for each_folder in all_folders_loci: # Looping through each of the locus folders
		outputfile_name = str(each_folder) + '_clustered.txt' # "each_folder" is also the name of the locus
		outputfile = open(outputfile_name, 'w')
		os.chdir(each_folder)
		bcodesForThisLocus = glob.glob("*")
		#print "glob thinks there are these barcode folders present: ", bcodesForThisLocus, "\n\n"
		for bcode_folder in bcodesForThisLocus: # have to go into each barcode folder in each locus folder
			if os.path.isdir(bcode_folder): # the glob might have found some files as well as folders	
				os.chdir(bcode_folder)
				files_to_add = glob.glob("*Ss3.fa")
				for file in files_to_add: 
					shutil.copyfileobj(open(file,'rb'), outputfile) #Add each file to the final output		
				os.chdir('..')
		os.chdir('..')
		outputfile.close()


	return LocusTaxonCountDict_clustd

ppp_location = os.path.dirname(os.path.abspath( __file__ ))
Usearch = ppp_location + '/' + 'Dependencies/usearch8.1.1756'

if len(sys.argv) < 6:
	sys.exit("""

Use this script to split an annotated fasta file based on taxon, barcode, loci, or group. 

Usage: ./purc_recluster.py annoated_file clustID1 clustID2 clustID sizeThreshold1 sizeThreshold2
Example: ./purc_recluster.py purc_run_3_annotated.fa 0.997 0.995 0.99 1 4

Note: 
(1) clustID1-3 : The similarity criterion for the first, second and third USEARCH clustering
(2) sizeThreshold : The min. number of sequences/cluster necessary for that cluster to be retained (set to 2 to remove singletons, 3 to remove singletons and doubles, etc)

	""")

log = open('purc_log.txt', 'w')
annotated_file = sys.argv[1]
clustID = float(sys.argv[2])
clustID2 = float(sys.argv[3])
clustID3 = float(sys.argv[4])
sizeThreshold = int(sys.argv[5])
sizeThreshold2 = int(sys.argv[6])

LocusTaxonCountDict_clustd = ClusterDechimera(annotated_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2)



taxon_list = []
locus_list = []
for taxon_locus in LocusTaxonCountDict_clustd.keys():
	taxon_list.append(taxon_locus[0])
	locus_list.append(taxon_locus[1])
taxon_list = set(taxon_list)
locus_list = set(locus_list)



count_output = open('purc_cluster_counts.xls', 'w')
count_output.write('\n**Final clustered sequences per accession per locus**\n')
log.write('\n**Final clustered sequences per accession per locus**\n')

count_output.write('\t' + '\t'.join(locus_list) + '\n')	
log.write('\t' + '\t'.join(locus_list) + '\n')	

for each_taxon in set(taxon_list):
	#print each_taxon, '\t',
	count_output.write(each_taxon + '\t')		
	log.write(each_taxon + '\t')		

	for each_locus in locus_list:
		try:
			#print LocusTaxonCountDict_clustd[each_taxon, each_locus], '\t',
			count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
			log.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
		except:
			#print '0', '\t',
			count_output.write('0\t') 
			log.write('0\t') 
	#print
	count_output.write('\n')		
	log.write('\n')		

#print 
#print '\n**Allele/copy/cluster/whatever count by locus**'
count_output.write('\n**Allele/copy/cluster/whatever count by locus**\n')	
log.write('\n**Allele/copy/cluster/whatever count by locus**\n')	

# I think this was breaking if a locus had no sequences, and thus that file is not created. Going to try "try"
for each_locus in locus_list:
	file_name = str(each_locus) + '_clustered.txt'
	try: #I'm hoping that this will stop the program from crashing if a locus has no sequences
		seq_no = len(parse_fasta(file_name))
		#print '\t', each_locus, ':', seq_no
		count_output.write(str(each_locus) + '\t' + str(seq_no) + '\n')
		log.write(str(each_locus) + '\t' + str(seq_no) + '\n')
	except:
		#print '\t', each_locus, ':', 0
		count_output.write(str(each_locus) + '\t0\n')			
		log.write(str(each_locus) + '\t0\n')			












