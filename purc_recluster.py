#!/usr/bin/env python

logo = """
-------------------------------------------------------------
|                       purc_recluster                      |
|        Pipeline for Untangling Reticulate Complexes       |
|                        version 1.0                        | 
|            https://bitbucket.org/crothfels/ppp            |
|															|
|                 Fay-Wei Li & Carl J Rothfels              |
|           see purc.py script for more information         |
-------------------------------------------------------------
""" 

usage = """

Use this script to recluster the alleles/homeologs from a previous PURC run. 

Example: ./purc_recluster.py -f purc_3_annotated.fa -o recluster -c 0.997 0.995 0.99 0.997 -s 1 4 --clean

"""

citation = """
This script relies heavily on USEARCH, MUSCLE, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of USEARCH, etc., used).

PURC: 
-Rothfels, C.J., K.M. Pryer, and F-W. Li. 2016. Next-generation polyploid 
phylogenetics: rapid resolution of hybrid polyploid complexes using PacBio 
single-molecule sequencing. New Phytologist

USEARCH/UCLUST: 
-Edgar, R.C. 2010. Search and clustering orders of magnitude faster than BLAST. 
Bioinformatics 26(19), 2460-2461.

UCHIME:
-Edgar, R.C., B.J. Haas, J.C. Clemente, C. Quince, R. Knight. 2011. 
UCHIME improves sensitivity and speed of chimera detection, Bioinformatics 27(16), 2194-2200.

Cutadapt:
-Martin, M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. 
EMBnet.journal 17:10-12.

MUSCLE:
-Edgar, R.C. 2004. MUSCLE: Multiple sequence alignment with high accuracy and high throughput. 
Nucleic Acids Research 32:1792-1797.

BLAST: 
-Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, et al. 2009. 
BLAST+: Architecture and applications. BMC Bioinformatics 10: 421.
"""


import sys
import os
import re
import subprocess
import glob
import shutil
import time
import datetime
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo

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
				# The locus names are the same as locus_folder
			except:
				LocusTaxonCountDict_unclustd[str(each_seq.id).split('|')[0], str(each_seq.id).split('|')[1]] = 1		
		
		os.chdir('..')
	if split_by == "locus":
		return splits_count_dic, LocusTaxonCountDict_unclustd #as {'BC01': 150, 'BC02': 156} for example
	else:
		return splits_count_dic

def clusterIt(file, clustID, round, previousClusterToCentroid_dict, verbose_level=1):
	"""The clustering step, using the clustID value"""
	outFile = re.sub(r"(.*)\.fa", r"\1C%s_%s.fa" %(round, clustID), file) # The rs indicate "raw" and thus python's escaping gets turned off
	outClustFile = re.sub(r"(.*).fa", r"\1clusts%s.uc" %(round), file)
	if round == 1:
		usearch_cline = "%s -cluster_fast %s -id %f -gapopen 3I/1E -consout %s -uc %s -sizeout" % (Usearch, file, clustID, outFile, outClustFile) 
		#usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -sortedby other -centroids %s -uc %s -sizeout" % (Usearch, file, clustID, outFile, outClustFile) 
		#usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -usersort -consout %s -uc %s -sizeout" % (Usearch, file, clustID, outFile, outClustFile) # Usearch 7   
        # Can add in "-cons_truncate" to the usearch call, if the primer removal isn't effective, but note some problems with partial sequences results.
	elif round > 1:
		usearch_cline = "%s -cluster_fast %s -id %f -gapopen 3I/1E -consout %s -uc %s -sizein -sizeout" % (Usearch, file, clustID, outFile, outClustFile)
		#usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -sortedby other -centroids %s -uc %s -sizein -sizeout" % (Usearch, file, clustID, outFile, outClustFile)
		#usearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -usersort -consout %s -uc %s -sizein -sizeout" % (Usearch, file, clustID, outFile, outClustFile) # Usearch 7	
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		#print '\n**Usearch-clustering output on', file, '**\n'
		#print err
		log.write('\n**Usearch-clustering output on' + str(file) + '**\n')
		log.write(str(err))
	
	uc = open(outClustFile, 'rU')
	#if round == 1: 
	ClusterToCentroid_dict = {}
	#	global ClusterToCentroid_dict
	for line in uc:
		if line.startswith('C'):
			cluster_name = 'Cluster' + str(line.split('\t')[1])
			centroid_seq_name = line.split('\t')[-2]
			
			if round > 1:
				ClusterToCentroid_dict[cluster_name] = previousClusterToCentroid_dict[centroid_seq_name.split(';')[0]]
			elif round == 1:
				ClusterToCentroid_dict[cluster_name] = centroid_seq_name

	if round == 4:
		try:
			clustered_seqs = parse_fasta(outFile)
			renamed_clustered_seqs = open('temp', 'a')
			for seq in clustered_seqs: # seq as 'Cluster15;size=1;'
				#new_seq_id = ClusterToCentroid_dict[str(seq.id).split(';')[0]]
				new_seq_id = str(seq.id).replace(str(seq.id).split(';')[0], ClusterToCentroid_dict[str(seq.id).split(';')[0]])
				#new_seq_id = ClusterToCentroid_dict[str(seq.id).split(';')[0]] + ';' + str(seq.id).split(';')[0:]
				#print str(seq.id).split(';')
				renamed_clustered_seqs.write('>' + new_seq_id + '\n' + str(seq.seq) + '\n')
			os.remove(outFile)
			os.rename('temp', outFile)
		except:
			log.write(outFile + ' is empty; perhaps sizeThreshold too high\n')

	return outFile, ClusterToCentroid_dict, outClustFile

# This function looks for PCR chimeras -- those formed within a single amplicon pool
def deChimeIt(file, round, abskew=1.9, verbose_level=0):
	"""Chimera remover"""
	outFile = re.sub(r"(.*)\.fa", r"\1dCh%s.fa" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
	outFile_uchime = re.sub(r"(.*)\.fa", r"\1dCh%s.uchime" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
	usearch_cline = "%s -uchime_denovo %s -abskew %s -nonchimeras %s -uchimeout %s" % (Usearch, file, abskew, outFile, outFile_uchime)
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level in [1, 2]:
		log.write('\n**Uchime output on' + str(file) + '**\n')
		log.write(str(err))

	# count number of chimera seq found by parsing 'outFile_uchime'
	chimera_count = 0
	uchime_out = open(outFile_uchime, 'rU')
	for line in uchime_out:
		line = line.strip('\n')
		if line.split('\t')[-1] == 'Y':
			chimera_count = chimera_count + 1
	return outFile, chimera_count

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
	
def sortIt_size(file, thresh, round, verbose_level=0):
	"""Sorts clusters by size, and removes those that are smaller than a particular size
	(sent as thresh -- ie, sizeThreshold). 
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

def align_and_consensus(inputfile, output_prefix):
	output_alignment = output_prefix.split(';')[0] + '_aligned.fa'
	output_consensus = output_prefix.split(';')[0] + '_consensus.fa'
	muscle_cline = '%s -in %s -out %s' % (Muscle, inputfile, output_alignment)
	process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 

	alignment = AlignIO.read(output_alignment, 'fasta')
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.gap_consensus(ambiguous='N',threshold=0.51)
	output = open(output_consensus, 'w')
	output.write('>' + output_prefix + '\n' + str(consensus).replace('-','') + '\n')
	return

def IterativeClusterDechimera(annotd_seqs_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2, verbose_level = 1): # M_p_barcode was set to FALSE, whcih led to splitting by taxon instead of barcode
	## Split sequences into separate files/folders for each locus ##
	sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
	locusCounts, LocusTaxonCountDict_unclustd = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus") 

	sys.stderr.write('Clustering/dechimera-izing seqs...\n')
	log.write('#Sequence clustering/dechimera-izing#\n')
	all_folders_loci = locusCounts.keys() # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
		# folders) and they correspond to the counts for each.
	LocusTaxonCountDict_clustd = {}
	LocusTaxonCountDict_chimera = {}
	
	## Go through each locus ##
	for locus_folder in all_folders_loci: # locus_folder = locus name
		os.chdir(locus_folder)
		sys.stderr.write('\nWorking on: ' + locus_folder + '...\n')
		if verbose_level in [1,2]:
			log.write('\nWorking on ' + str(locus_folder) + ' ...\n')

		if not os.stat(locus_folder + ".fa").st_size == 0: # ie, the file is not empty
			## Split sequences into separate taxon folders ##
			taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon")					
			all_folders_taxon = taxonCounts.keys()

			for taxon_folder in all_folders_taxon:
				if verbose_level in [1,2]:
					log.write("Working on " + taxon_folder + '\n')
				os.chdir(taxon_folder)
				
				if verbose_level in [1,2]:
					log.write("\tFirst clustering\n")
					log.write("\nAttempting to sort: " + taxon_folder + ".fa\n")

				sorted_length = sortIt_length(file = taxon_folder + ".fa", verbose_level = verbose_level)
				clustered1, previousClusterToCentroid_dict, outClustFile1 = clusterIt(file = sorted_length, previousClusterToCentroid_dict = '', clustID = clustID, round = 1, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tFirst chimera slaying expedition\n")
				deChimered1, chimera_count1 = deChimeIt(file = clustered1, round = 1, abskew = abskew, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tSecond clustering\n")
				sorted_size1 = sortIt_size(file = deChimered1, thresh = sizeThreshold, round = 1, verbose_level = verbose_level)
				clustered2, previousClusterToCentroid_dict, outClustFile2 = clusterIt(file = sorted_size1, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID2, round = 2, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tSecond chimera slaying expedition\n")
				deChimered2, chimera_count2 = deChimeIt(file = clustered2, round = 2, abskew = abskew, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tThird clustering\n")
				sorted_size2 = sortIt_size(file = deChimered2, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
				clustered3, previousClusterToCentroid_dict, outClustFile3 = clusterIt(file = sorted_size2, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID3, round = 3, verbose_level = verbose_level)
				
				if verbose_level in [1,2]:
					log.write("\tThird chimera slaying expedition\n")
				deChimered3, chimera_count3 = deChimeIt(file = clustered3, round = 3, abskew = abskew, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\Forth clustering\n")
				sorted_size3 = sortIt_size(file = deChimered3, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
				clustered4, previousClusterToCentroid_dict, outClustFile4 = clusterIt(file = sorted_size3, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID3, round = 4, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tThird chimera slaying expedition\n")
				deChimered4, chimera_count4 = deChimeIt(file = clustered4, round = 4, abskew = abskew, verbose_level = verbose_level)

				sorted_size4 = sortIt_size(file = deChimered4, thresh = sizeThreshold2, round = 4, verbose_level = verbose_level)

				try:
					clustered_seq_file = parse_fasta(sorted_size4)
					for each_seq in clustered_seq_file:
						#taxon_name = str(each_seq.id).split('|')[0].split('=')[-1] # for example, get C_dia_5316 from centroid=centroid=C_dia_5316|ApP|C|BC02|_p0/158510/ccs;ee=1.9;;seqs=6;seqs=18;size=27;					
						try:
							LocusTaxonCountDict_clustd[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
							# The locus names are the same as locus_folder
						except:
							LocusTaxonCountDict_clustd[taxon_folder, locus_folder] = 1		
				except:
					if verbose_level in [1,2]:
						log.write(str(sorted_size4) + 'is an empty file\n')

				### Collect all sequences from each cluster and re-consensus ###
				# Go through the first clustering uc file
				ClusterToSeq_dict1 = {}
				for line in open(outClustFile1, 'rU'):	
					line = line.strip('\n')
					if line.startswith('H') or line.startswith('C'):
						key = 'Cluster' + line.split('\t')[1]
						seq = line.split('\t')[8]
						try:
							ClusterToSeq_dict1[key].append(seq)
						except:
							ClusterToSeq_dict1[key] = [seq]
				# Go through the second clustering uc file
				ClusterToSeq_dict2 = {}
				for line in open(outClustFile2, 'rU'):	
					line = line.strip('\n')
					if line.startswith('H') or line.startswith('C'):
						key = 'Cluster' + line.split('\t')[1]
						seqs = ClusterToSeq_dict1[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
						for seq in seqs:
							try:
								ClusterToSeq_dict2[key].append(seq)
							except:
								ClusterToSeq_dict2[key] = [seq]
				# Go through the third clustering uc file
				ClusterToSeq_dict3 = {}
				for line in open(outClustFile3, 'rU'):	
					line = line.strip('\n')
					if line.startswith('H') or line.startswith('C'):
						key = 'Cluster' + line.split('\t')[1]
						seqs = ClusterToSeq_dict2[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
						for seq in seqs:
							try:
								ClusterToSeq_dict3[key].append(seq)
							except:
								ClusterToSeq_dict3[key] = [seq]
				# Go through the forth clustering uc file
				ClusterToSeq_dict4 = {}
				for line in open(outClustFile4, 'rU'):	
					line = line.strip('\n')
					if line.startswith('H') or line.startswith('C'):
						key = 'Cluster' + line.split('\t')[1]
						seqs = ClusterToSeq_dict3[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
						for seq in seqs:
							try:
								ClusterToSeq_dict4[key].append(seq)
							except:
								ClusterToSeq_dict4[key] = [seq]
				
				## Align and re-consensus of all constituent seq for each cluster ##
				SeqDict = SeqIO.index(taxon_folder + ".fa", 'fasta')
				for each_cluster in ClusterToSeq_dict4:
					if len(ClusterToSeq_dict4[each_cluster]) >= int(sizeThreshold2):
						cluster_seq_file = open(each_cluster, 'w')
						for seq in ClusterToSeq_dict4[each_cluster]:
							cluster_seq_file.write('>' + seq + '\n' + str(SeqDict[seq].seq) + '\n')
						cluster_seq_file.close()
						new_seq_name = taxon_folder + '_' + each_cluster + ';size=' + str(len(ClusterToSeq_dict4[each_cluster])) + ';'
						align_and_consensus(each_cluster, new_seq_name)
				
				## Put all consensus seq into one file *_Cluster_Finalconsensus.fa ##
				all_consensus_seq = open(taxon_folder + '_Cluster_Finalconsensus.fa', 'w')
				files_to_add_reconsensus = glob.glob("*_consensus.fa")
				for file in files_to_add_reconsensus: 
					shutil.copyfileobj(open(file,'rb'), all_consensus_seq) #Add each file to the final output
				all_consensus_seq.close()

				## Do final clustering and chimera-killing ##
				usearch_cline = "%s -sortbysize %s -fastaout %s" %(Usearch, taxon_folder + '_Cluster_Finalconsensus.fa', taxon_folder + '_Cluster_FinalconsensusSs.fa')
				process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
				(out, err) = process.communicate() #the stdout and stderr
				
				usearch_cline = "%s -cluster_fast %s -id %f -gapopen 3I/1E -consout %s -uc %s -sizein -sizeout" % (Usearch, taxon_folder + '_Cluster_FinalconsensusSs.fa', clustID4, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.uc')
				process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
				(out, err) = process.communicate() #the stdout and stderr

				usearch_cline = "%s -uchime_denovo %s -abskew %s -nonchimeras %s -uchimeout %s" % (Usearch, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.fa', abskew, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.uchime')
				process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
				(out, err) = process.communicate() #the stdout and stderr
				
				## Count number of chimera seq found by parsing 'outFile_uchime' ##
				chimera_count5 = 0
				uchime_out = open(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.uchime', 'rU')
				for line in uchime_out:
					line = line.strip('\n')
					if line.split('\t')[-1] == 'Y':
						chimera_count5 = chimera_count5 + 1
				LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = [chimera_count1, chimera_count2, chimera_count3, chimera_count4, chimera_count5]

				## Rename sequences in the final fasta: add taxon name ##
				sed_cmd = "sed 's/>/>%s_/g' %s > %s" % (taxon_folder, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh_renamed.fa')
				process = subprocess.Popen(sed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
				(out, err) = process.communicate() #the stdout and stderr
				
				## Remove intermediate files if requsted ##
				if args.clean:
					filt_to_remove = list(set(glob.glob('*')) - set(glob.glob('*_renamed.fa')) - set(glob.glob(taxon_folder+'.fa')))
					for file in filt_to_remove:
						os.remove(file)

				os.chdir("..") # To get out of the current barcode folder and ready for the next one
		os.chdir("..") # To get out of the current locus folder and ready for the next one
	log.write('\t...done\n\n')	

	## Put all the sequences together ##
	sys.stderr.write('\rPutting all the sequences together......\n\n')
	for locus_folder in all_folders_loci: # Looping through each of the locus folders
		#outputfile_name = str(locus_folder) + '_clustered.txt' # "locus_folder" is also the name of the locus
		#outputfile = open(outputfile_name, 'w')
		outputfile_name_reconsensus = str(locus_folder) + '_clustered_reconsensus.fa' 
		outputfile_reconsensus = open(outputfile_name_reconsensus, 'w')

		os.chdir(locus_folder)
		taxonForThisLocus = glob.glob("*")
		for taxon_folder in taxonForThisLocus: # have to go into each barcode folder in each locus folder
			if os.path.isdir(taxon_folder): # the glob might have found some files as well as folders	
				os.chdir(taxon_folder)
				#files_to_add = glob.glob("*Ss4.fa")				
				#for file in files_to_add: 
				#	shutil.copyfileobj(open(file,'rb'), outputfile) #Add each file to the final output
				files_to_add_reconsensus = glob.glob("*_Cluster_FinalconsensusSsC*_renamed.fa")
				for file in files_to_add_reconsensus: 
					shutil.copyfileobj(open(file,'rb'), outputfile_reconsensus) #Add each file to the final output

				os.chdir('..')
		os.chdir('..')
		outputfile_reconsensus.close()

	return LocusTaxonCountDict_clustd, LocusTaxonCountDict_unclustd, LocusTaxonCountDict_chimera

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

# Parse the arguments #
parser = argparse.ArgumentParser(description=usage)

required = parser.add_argument_group('required arguments')
required.add_argument('-f','--annotated_file', type=str, action="store", metavar='\b', required=True,
               help='The annotated fasta file, generated by purc')
required.add_argument('-o','--output_folder', type=str, action="store", metavar='\b', required=True,
               help='The output folder')
required.add_argument('-c','--clustering_identities', type=float, nargs=4, metavar='\b', required=True,
               help='The similarity criterion for the first, second, third and forth USEARCH clustering',
               default=[0.997, 0.995, 0.99, 0.997])
required.add_argument('-s','--size_threshold', type=int, nargs=2, metavar='\b', required=True,
                    help='The min. number of sequences/cluster necessary for that cluster to be retained (set to 2 to remove singletons, 3 to remove singletons and doubles, etc)',
                    default=[1, 4])

parser.add_argument('-a','--abuncdance_skew', type=float, nargs=1, metavar='\b',
                    help='The parameter to control chimera-killing; the default is 1.9',
                    default=1.9)
parser.add_argument('--clean', action="store_true",
                    help='Remove the intermediate files')

args = parser.parse_args()

annotated_file = args.annotated_file
masterFolder = args.output_folder
clustID = args.clustering_identities[0]
clustID2 = args.clustering_identities[1]
clustID3 = args.clustering_identities[2]
clustID4 = args.clustering_identities[3]
sizeThreshold = args.size_threshold[0]
sizeThreshold2 = args.size_threshold[1]
abskew = args.abuncdance_skew[0]

purc_location = os.path.dirname(os.path.abspath( __file__ ))
Usearch = purc_location + '/' + 'Dependencies/usearch8.1.1756'
Muscle = purc_location + '/' + 'Dependencies/muscle3.8.31'

# Check if dependencies are in place
sys.stderr.write('Checking dependencies...\n')
if not os.path.isfile(Usearch):
	sys.exit("Error: couldn't find the Usearch executable")
if not os.path.isfile(Muscle):
	sys.exit("Error: couldn't find the Muscle executable")

# Check if muscle can be executed
muscle_cline = '%s -version' % (Muscle)
process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
(out, err) = process.communicate() #the stdout and stderr
if not str(out).startswith('MUSCLE'):
	sys.exit("Error: could not execute Muscle")

# Check if Usearch can be executed
usearch_cline = '%s -version' % (Usearch)
process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
(out, err) = process.communicate() #the stdout and stderr
if not str(out).startswith('usearch'):
	sys.exit("Error: could not execute Usearch")

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
log.write('purc_recluster.py ' + annotated_file + ' ' + masterFolder + ' ' + str(clustID) + ' ' + str(clustID2) + ' ' + str(clustID3) + ' ' + str(sizeThreshold) + ' ' + str(sizeThreshold2) + '\n\n')
if not os.path.isfile(annotated_file):
	sys.exit("Error: could not find " + annotated_file)

## Recluster and redechimera ##
LocusTaxonCountDict_clustd, LocusTaxonCountDict_unclustd, LocusTaxonCountDict_chimera = IterativeClusterDechimera('../'+annotated_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2)

taxon_list = []
locus_list = []
for taxon_locus in LocusTaxonCountDict_clustd.keys():
	taxon_list.append(taxon_locus[0])
	locus_list.append(taxon_locus[1])
taxon_list = set(taxon_list)
locus_list = set(locus_list)

## Producing a summary ##
count_output = open('purc_cluster_counts.xls', 'w')

# Output read count per taxon per locus
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

# Output clustered seq count
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

# Output cluster count per locus 
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

# Output chimeric seq count
count_output.write('\n**Chimeric clusters/sequence count by locus**\n')	
for each_locus in locus_list:
	count_output.write(each_locus + '\n')		
	#log.write(each_locus + '\t')		
	for each_taxon in set(taxon_list):
		try:
			count_output.write('\t' + each_taxon + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][0]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][1]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][2]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][3]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][4]))		
		except:
			count_output.write('\t' + each_taxon)	
		
		count_output.write('\n')		
	#log.write('\n')		

## Aligning the sequences ##
fastas = glob.glob("*_renamed.fa")
for file in fastas:
	sys.stderr.write("Aligning " + file + "\n")
	log.write("Aligning " + file + "\n")
	outFile = muscleIt(file)
	
fastas = glob.glob("*_clustered_reconsensus.fa")
for file in fastas:
	sys.stderr.write("Aligning " + file + "\n")
	log.write("Aligning " + file + "\n")
	outFile = muscleIt(file)









