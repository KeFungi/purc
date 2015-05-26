#!/usr/bin/env python

#PUSH PUSH to see PPP
#TODO - create a mode that just annotates (no splitting or clustering)

# TODO - Maybe the secret would be to include "dummy" entries in the maps. Ie., if there's no groupB for a given
# barcode, map groupB for that barcode to something like "MisMatch". Then they'll show up in the trees
# etc and can be dealt with.

# TODO re:mode 2. The naming of files and folders needs to get clean up. Currently this mode requires that the output folder
# be already created (this is ok, I suppose), and that the input seqence file be named output_prefix + annotated.fa? It's a
# little weirder than that -- when I had the Input_sequence_file set to PBR2_C_3_annotated.fa and the Output_prefix set to
# splitByBC, I got the following error: No such file or directory: 'splitByBC_3_annotated.fa'. When I renamed the input seq
# file to that, and stuck inside the splitByBC folder, it worked

# TODO re:mode 3. Currently mode 2 creates subdirectories for each of the splits it makes, but mode 3 wants all the files
# to be in a single directory. Do we want to fix this? Could leave it as is, with instructions..

# TODO: add in a cleanUpFiles option in the .conf such that if it is set it yes it erases all the S1C10.999CH1 etc files when they are no longer needed
# if we do this (perhaps worth doing regardless), we should get PUSSwhatever to spit out the settings used (specifically the clustering settings)
# so that that information doesn't get lost (the researcher can go back and figure out what they did...)

print
print """
-------------------------------------------------------------
|                             PPP                           |
|                 PacBio Polyploidy Pipeline                |
|                        version 1.90                       | 
|                                                           |
|                 Fay-Wei Li & Carl J Rothfels              |
-------------------------------------------------------------
"""
print 

usage = """
You need to provide a configuration file.

Usage: ./ppp.py configuration_file > out
Example: ./ppp.py ppp_configuration.txt > summary.txt
For more info, try: ./ppp.py -help

"""

citation = """
This script relies heavily on USEARCH, MUSCLE, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of USEARCH, etc., used).

PPP: 
-Awesome paper by carl and fay-wei. Awesome journal. Awesome page numbers.

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

import re
import sys
import os
import glob
import subprocess
import shutil
from Bio import SeqIO

def parse_fasta(infile):
	"""Reads in a fasta, returns a list of biopython seq objects"""
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq] 

def count_seq_from_fasta(infile):
	cmdLine = "grep '>' %s | wc -w" % infile
	process = subprocess.Popen(cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate()	
	out = out.strip('\n').replace(' ', '')
	return out

def ReverseComplement(seq):
	"""Returns reverse complement sequence, ignores gaps"""
	seq = seq.replace(' ','') # Remove spaces
	seq = seq[::-1] # Reverse the sequence
	basecomplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y':'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W', 'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B'} # Make a dictionary for complement
	letters = list(seq) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters) 

def makeBlastDB(inFileName, outDBname): #FWL do we need to include instructions for users to get the makeblastdb command to work?
	"""makes a blast database from the input file"""
	makeblastdb_cmd = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (inFileName, outDBname)
	process = subprocess.Popen(makeblastdb_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate()
	print out, err
	return

def BlastSeq(inputfile, outputfile, databasefile, evalue=0.0000001, max_target=1, outfmt='6 qacc sacc nident mismatch length pident bitscore'):	
	"""Calls blastn, the output format can be changed by outfmt. Requires the blast database to be made already"""
	blastn_cLine = "blastn -query %s -task blastn -num_threads 8 -db %s -out %s -evalue %s -max_target_seqs %d -outfmt '%s'" % (inputfile, databasefile, outputfile, evalue, max_target, outfmt)
	os.popen(blastn_cLine)	
	return

def CheckChimericLoci(inputfile_raw_sequences, outputfile_blast, outputfile_goodSeqs, outputfile_chimeras, databasefile, SeqDict, SplitChimera=False):
	"""Blast each input sequences to the reference database to detect chimeric sequences (i.e. a sequence hit to two different loci) 
	Return "chimera_dict", in which the sequence name is the key and [locus_name1, locus_name2] is the value
	If SplitChimera = True, then will split the chimeric sequence into two loci, based on the coordinates returned from BLAST. 
		NOTE: to ensure the BC are split together with each locus, the orientation/strand also matters. 
	"""
	BlastSeq(inputfile_raw_sequences, outputfile_blast, databasefile, evalue=1e-100, max_target=100, outfmt='6 qacc sacc length pident evalue qstart qend qlen sstrand')
	
	chimera_blast = open(outputfile_blast, 'rU') # Read the blast result
	chimeras = open(outputfile_chimeras, 'w')
	non_chimeras = open(outputfile_goodSeqs, 'w')

	loci_info_dict = {}
	chimera_dict = {}
	
	for each_rec in chimera_blast:
		each_rec = each_rec.strip('\n')
		seq_name = each_rec.split('\t')[0]
		First = True
		if seq_name not in chimera_dict.keys():
			locus_name = each_rec.split('\t')[1].split('_')[0].upper() #upper to avoid "ApP" != "APP" issue
			locus_start = each_rec.split('\t')[5]
			locus_end = each_rec.split('\t')[6]
			locus_direction = each_rec.split('\t')[-1]
			if seq_name not in loci_info_dict.keys():
				loci_info_dict[seq_name] = [(locus_name, int(locus_start), int(locus_end), locus_direction)]

			else:
				if locus_name != loci_info_dict[seq_name][0][0]:
					chimera_dict[seq_name] = [locus_name, loci_info_dict[seq_name][0][0]]
					if First:
						loci_info_dict[seq_name].append((locus_name, int(locus_start), int(locus_end), locus_direction))
						First = False

	non_chimeras_list = list(set(list(SeqDict.keys())) - set(list(chimera_dict.keys())))
	for each_non_chimera in non_chimeras_list:
		non_chimeras.write('>' + each_non_chimera + '\n' + str(SeqDict[each_non_chimera].seq) + '\n')

	if SplitChimera:
		for each_chimera in chimera_dict:
			loci_info_list = sorted(loci_info_dict[each_chimera], key=lambda LIST: LIST[1]) # loci_info_list as [('IBR', 41, 906, 'plus'), ('APP', 956, 1920, 'minus')]

			new_seq_name1 = each_chimera + 'ERRchimera_' + loci_info_list[0][0] + 'of' + loci_info_list[0][0] + '+' + loci_info_list[1][0]
			new_seq_name2 = each_chimera + 'ERRchimera_' + loci_info_list[1][0] + 'of' + loci_info_list[0][0] + '+' + loci_info_list[1][0]
			
			if loci_info_list[0][3] == 'plus' and loci_info_list[1][3] == 'plus':
				locus_seq1 = str(SeqDict[each_chimera].seq)[:loci_info_list[0][2]]
				locus_seq2 = str(SeqDict[each_chimera].seq)[loci_info_list[0][2]:]

			elif loci_info_list[0][3] == 'plus' and loci_info_list[1][3] == 'minus':
				locus_seq1 = str(SeqDict[each_chimera].seq)[:loci_info_list[0][2]]
				locus_seq2 = str(SeqDict[each_chimera].seq)[loci_info_list[0][2]:]
			
			elif loci_info_list[0][3] == 'minus' and loci_info_list[1][3] == 'plus':
				midpoint = int((loci_info_list[0][2] + loci_info_list[1][1]) / 2)
				locus_seq1 = str(SeqDict[each_chimera].seq)[:midpoint]
				locus_seq2 = str(SeqDict[each_chimera].seq)[midpoint:]

			elif loci_info_list[0][3] == 'minus' and loci_info_list[1][3] == 'minus':
				locus_seq1 = str(SeqDict[each_chimera].seq)[:loci_info_list[1][1]]
				locus_seq2 = str(SeqDict[each_chimera].seq)[loci_info_list[1][1]:]

			#for each_locus in sorted(loci_info_dict[each_chimera], key=lambda list: list[1]): #each_locus as (locus_name, start_posi, end_posi)
				#print each_locus
				#new_seq_name = each_chimera + 'ERRchimera_' + each_locus[0] + 'of' + chimera_dict[each_chimera][0] + '+' + chimera_dict[each_chimera][1]
				#locus_seq = str(SeqDict[each_chimera].seq)[int(each_locus[1]):int(each_locus[2])]
			non_chimeras.write('>' + new_seq_name1 + '\n' + locus_seq1 + '\n')
			non_chimeras.write('>' + new_seq_name2 + '\n' + locus_seq2 + '\n')

	for each_chimera in chimera_dict:
		new_seq_name = each_chimera + 'ERRchimera_' + chimera_dict[each_chimera][0] + '+' + chimera_dict[each_chimera][1]
		chimeras.write('>' + new_seq_name + '\n' + str(SeqDict[each_chimera].seq) + '\n')


	return chimera_dict
	
def DeBarcoder(inputfile_raw_sequences, databasefile, SeqDict, Output_folder, Output_prefix):
	"""Blasts the raw sequences against the barcode blast database, identifies the barcode, adds the barcode ID to the 
	sequence name, removes the barcode from sequence; returns a list of barcode names"""
	
	BlastSeq(inputfile_raw_sequences, Output_folder + '/blast_barcode_out.txt', databasefile, evalue=1, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')
	
	bc_blast = open(Output_folder + '/blast_barcode_out.txt', 'rU') # Read the blast result
	bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'w') # For writing the de-barcoded sequences
	bc_leftover = open(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'w') # For saving those without barcodes
	bc_toomany = open(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa', 'w') # For saving those more than one barcode

	seq_withbc_list = [] # A list containing all the seq names that have barcodes
	seq_withbc_morethanone_list = [] # A list containing all the seq names that have more than one barcode
	seq_withoutbc_list = [] # A list containing all the seq names that do not have barcode identified by BLAST

	barcode_info_dict = {} # {seq_name1: [BC01, 0, 12], seq_name2: [BC08, 0, 12]}; barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi]
	
	# Go through the blast output file, and complete the barcode_info_dict, seq_withbc_list, and seq_withbc_morethanone_list
	for each_rec in bc_blast:
		each_rec = each_rec.strip('\n')
		seq_name = each_rec.split('\t')[0]
		barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
		barcode_start_posi = int(each_rec.split('\t')[5])
		barcode_end_posi = int(each_rec.split('\t')[6])		
		
		if seq_name not in barcode_info_dict.keys():
			barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi]
			seq_withbc_list.append(seq_name)
		else: # means that this seq has more than one barcode, then take out this seq record from seq_withbc_list, but append it to seq_withbc_morethanone_list
			del barcode_info_dict[seq_name]
			seq_withbc_list.remove(seq_name)
			seq_withbc_morethanone_list.append(seq_name)

	# De-barcode and write sequences
	for each_seq in seq_withbc_list:
		new_seq_name = str(barcode_info_dict[each_seq][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name

		#check the orientation of the sequence; if the barcode is in the 3' end, reverse complement the seq
		if barcode_info_dict[each_seq][1] < 15: # The start position _should_ be 1, but this allows for some slop
			new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][2]:]) # "barcode_end_posi" to the end of the sequence
		elif barcode_info_dict[each_seq][1] > len(str(SeqDict[each_seq].seq))-40: # Those barcodes that are at the end of the sequences, so need to be reversecomplemented
			new_seq_trimmed = ReverseComplement(str(SeqDict[each_seq].seq[:barcode_info_dict[each_seq][1]-1])) # the beginning of the sequence to "barcode_start_posi" - 1
		else: # Those barcodes that are at the middle of the sequences
			new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][2]:])
			new_seq_name = new_seq_name + "ERRmidBC"
		bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
	
	# Write the sequences with multiple barcodes
	for each_seq in seq_withbc_morethanone_list:
		bc_toomany.write('>' + str(each_seq) + '\n' + str(SeqDict[each_seq].seq) + '\n')

	# Write the sequences without identified barcode to bc_leftover
	seq_withoutbc_list = list(set(list(SeqDict.keys())) - set(seq_withbc_list) - set(seq_withbc_morethanone_list))
	for seq_withoutbc in seq_withoutbc_list:
		bc_leftover.write('>' + str(seq_withoutbc) + '\n' + str(SeqDict[seq_withoutbc].seq) + '\n')
        
        
	bc_blast.close()
	bc_toomany.close()
	bc_leftover.close()
	bc_trimmed.close() #this is the file that now has all the sequences, labelled with the barcode, and the barcodes themselves removed
	
	return 

# This function looks for barcodes at the ends of the sequence, only. So it avoids internal barcodes (helpful in part because some of our loci had regions that were similar to some primer sequences)
def DeBarcoder_ends(SeqDict, databasefile, Output_folder, Output_prefix, search_range=25):
	#blast_outfile = open(outputfile, 'w')
	F_ends = open('tempF', 'w')
	R_ends = open('tempR', 'w')
	bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'w') # For writing the de-barcoded sequences
	bc_leftover = open(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'w') # For saving those without barcodes
	bc_toomany = open(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa', 'w') # For saving those more than one barcode

	for each_rec in sorted(SeqDict):
		#blastn_cLine = "blastn -db %s -task blastn -max_target_seqs 1 -outfmt '6 qacc sacc length pident evalue qstart qend qlen'" % databasefile
		#blastn_F = subprocess.Popen(blastn_cLine, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		#results_F, err = blastn_F.communicate(seq_to_search_F)		
		#blast_outfile.write(results_F)

		seq_to_search_F = str(SeqDict[each_rec].seq)[:search_range]
		F_ends.write('>' + str(each_rec) + '\n' + seq_to_search_F + '\n')

		seq_to_search_R = ReverseComplement(str(SeqDict[each_rec].seq))[:search_range]
		R_ends.write('>' + str(each_rec) + '\n' + seq_to_search_R + '\n')

	F_ends.close()
	R_ends.close()
	BlastSeq('tempF', Output_folder + '/blast_barcodeF_out.txt', databasefile, evalue=1, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')
	BlastSeq('tempR', Output_folder + '/blast_barcodeR_out.txt', databasefile, evalue=1, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')

	seq_withbc_list = [] # A list containing all the seq names that have barcodes
	seq_withbc_morethanone_list = [] # A list containing all the seq names that have more than one barcode
	seq_withoutbc_list = [] # A list containing all the seq names that do not have barcode identified by BLAST
	barcode_info_dict = {} # {seq_name1: [BC01, 0, 12], seq_name2: [BC08, 0, 12]}; barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi]

	bc_blast_F = open(Output_folder + '/blast_barcodeF_out.txt', 'rU')
	for each_rec in bc_blast_F:
		each_rec = each_rec.strip('\n')
		seq_name = each_rec.split('\t')[0]
		barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
		barcode_start_posi = int(each_rec.split('\t')[5])
		barcode_end_posi = int(each_rec.split('\t')[6])		
		if seq_name not in barcode_info_dict.keys():
			barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi, '+']
			seq_withbc_list.append(seq_name)
		else: # means that this seq has more than one barcode, then take out this seq record from seq_withbc_list, but append it to seq_withbc_morethanone_list
			del barcode_info_dict[seq_name]
			seq_withbc_list.remove(seq_name)
			seq_withbc_morethanone_list.append(seq_name)

	bc_blast_R = open(Output_folder + '/blast_barcodeR_out.txt', 'rU')
	for each_rec in bc_blast_R:
		each_rec = each_rec.strip('\n')
		seq_name = each_rec.split('\t')[0]
		barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
		barcode_start_posi = int(each_rec.split('\t')[5])
		barcode_end_posi = int(each_rec.split('\t')[6])		
		if seq_name not in barcode_info_dict.keys():
			barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi, '-']
			seq_withbc_list.append(seq_name)
		else: # means that this seq has more than one barcode, then take out this seq record from seq_withbc_list, but append it to seq_withbc_morethanone_list
			del barcode_info_dict[seq_name]
			seq_withbc_list.remove(seq_name)
			seq_withbc_morethanone_list.append(seq_name)

	# De-barcode and write sequences
	for each_seq in seq_withbc_list:
		new_seq_name = str(barcode_info_dict[each_seq][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name

		#check the orientation of the sequence; if the barcode is in the 3' end, reverse complement the seq
		if barcode_info_dict[each_seq][-1] == '+': 
			new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][2]:]) 
		elif barcode_info_dict[each_seq][-1] == '-': 
			new_seq_trimmed = ReverseComplement(str(SeqDict[each_seq].seq))[barcode_info_dict[each_seq][2]:]

		bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
	
	# Write the sequences with multiple barcodes
	for each_seq in seq_withbc_morethanone_list:
		bc_toomany.write('>' + str(each_seq) + '\n' + str(SeqDict[each_seq].seq) + '\n')

	# Write the sequences without identified barcode to bc_leftover
	seq_withoutbc_list = list(set(list(SeqDict.keys())) - set(seq_withbc_list) - set(seq_withbc_morethanone_list))
	for seq_withoutbc in seq_withoutbc_list:
		bc_leftover.write('>' + str(seq_withoutbc) + '\n' + str(SeqDict[seq_withoutbc].seq) + '\n')

	bc_blast_F.close()
	bc_blast_R.close()
	bc_toomany.close()
	bc_leftover.close()
	bc_trimmed.close() #this is the file that now has all the sequences, labelled with the barcode, and the barcodes themselves removed


def DeBarcoder_dual(inputfile_raw_sequences, databasefile, SeqDict):
	"""Blasts the raw sequences against the barcode blast database, identifies the barcode, adds the barcode ID to the 
	sequence name, removes the barcode from sequence; deal with barcodes at both primers"""
	BlastSeq(inputfile_raw_sequences, Output_folder + '/blast_barcode_out.txt', databasefile, evalue=1, max_target=2, outfmt='6 qacc sacc length pident bitscore qstart qend')
	
	bc_blast = open(Output_folder + '/blast_barcode_out.txt', 'rU') # Read the blast result
	bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'w') # For writing the de-barcoded sequences
	bc_leftover = open(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'w') # For saving those without barcodes
	bc_onlyF = open(Output_folder + '/' + Output_prefix + '_1_trashBin_onlyF_bc.fa', 'w')
	bc_onlyR = open(Output_folder + '/' + Output_prefix + '_1_trashBin_onlyR_bc.fa', 'w')
	bc_toomany = open(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa', 'w') # For saving those more than one barcode
	bc_invalid = open(Output_folder + '/' + Output_prefix + '_1_trashBin_invalid_bc.fa', 'w') # For saving those having FF or RR barcodes

	seq_withbc_list = [] # A list containing all the seq names that have barcodes
	seq_withoutbc_list = [] # A list containing all the seq names that do not have barcode identified by BLAST
	seq_onlyF_list = []
	seq_onlyR_list = []
	barcode_info_dict = {} # {seq_name1: [(BCF01, 0, 12), (BCR02, 459, 511)], seq_name2: [(BC08, 0, 12)]}; barcode_info_dict[seq_name] = [(barcodeF_name, barcode_start_posi, barcode_end_posi), (barcodeR_name, barcode_start_posi, barcode_end_posi)]
							# each seq has a list of tuples, each tuple contains the (barcode name, start, end)
	# Go through the blast output file, and complete the barcode_info_dict, and seq_withbc_list
	for each_rec in bc_blast:
		each_rec = each_rec.strip('\n')
		seq_name = each_rec.split('\t')[0]
		barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
		barcode_start_posi = int(each_rec.split('\t')[5])
		barcode_end_posi = int(each_rec.split('\t')[6])		
		seq_withbc_list.append(seq_name)
		# {seq_name1: [(BCF01, 0, 12), (BCR02, 459, 511)], seq_name2: [(BC08, 0, 12)]}
		try:
			barcode_info_dict[seq_name].append((barcode_name, barcode_start_posi, barcode_end_posi))
		except:
			barcode_info_dict[seq_name] = [(barcode_name, barcode_start_posi, barcode_end_posi)]
	
	# De-barcode and write sequences
	for each_seq in set(seq_withbc_list):
		# When more than three barcodes are present in seq
		if len(barcode_info_dict[each_seq]) >= 3:
			bc_toomany.write('>' + str(each_seq) + '\n' + str(SeqDict[each_seq].seq) + '\n')
			continue
		# When only one barcode is present in seq
		elif len(barcode_info_dict[each_seq]) == 1:	
			new_seq_name = str(barcode_info_dict[each_seq][0][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name
			# Check the orientation of the sequence; if the barcode is in the 3' end, reverse complement the seq		
			if barcode_info_dict[each_seq][0][1] < 5: # The start position _should_ be 1, but this allows for some slop
				new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][0][2]:]) # "barcode_end_posi" to the end of the sequence
			elif barcode_info_dict[each_seq][0][1] > len(str(SeqDict[each_seq].seq))-30: # Those barcodes that are at the end of the sequences, so need to be reversecomplemented
				new_seq_trimmed = ReverseComplement(str(SeqDict[each_seq].seq[:barcode_info_dict[each_seq][0][1]-1])) # the beginning of the sequence to "barcode_start_posi" - 1
			else: # Those barcodes that are at the middle of the sequences
				new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][0][2]:])
				new_seq_name = new_seq_name + "ERRmidBC"
			# Check where barcode is located, on F or R primer
			if barcode_info_dict[each_seq][0][0][2] == 'F':
				bc_onlyF.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
			elif barcode_info_dict[each_seq][0][0][2] == 'R':
				bc_onlyR.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')			
		# When both barcodes are present in seq
		elif len(barcode_info_dict[each_seq]) == 2:		
			# Check if barcodes are in the same category (F, R); don't want these 
			if barcode_info_dict[each_seq][0][0][2] == barcode_info_dict[each_seq][1][0][2]:
				#print 'Error: Invalid Dual-barcodes for', each_seq
				new_seq_name = str(barcode_info_dict[each_seq][0][0]) + '^' + str(barcode_info_dict[each_seq][1][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name
				bc_invalid.write('>' + new_seq_name + '\n' + str(SeqDict[each_seq].seq) + '\n')
				continue 

			# Re-order the list, so that the F barcode tuple is at the first position in the list
			if barcode_info_dict[each_seq][0][0][2] == 'R':
				barcode_info_dict[each_seq] = barcode_info_dict[each_seq][::-1]
			
			new_seq_name = str(barcode_info_dict[each_seq][0][0]) + '^' + str(barcode_info_dict[each_seq][1][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BCF01|BCR02|sequence_name

			# Trim the barcodes 
			# BCF - F primer - Seq - R primer - BCR
			if barcode_info_dict[each_seq][0][1] < 5: 
				new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][0][2]:barcode_info_dict[each_seq][1][1]-1]) 
				bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
			# BCR - R primer - Seq - F primer - BCF
			elif barcode_info_dict[each_seq][0][1] > len(str(SeqDict[each_seq].seq))-30: # When the seq is reverse complemented
				new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][1][2]:barcode_info_dict[each_seq][0][1]-1]) 
				new_seq_trimmed = ReverseComplement(new_seq_trimmed)		
				bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
			else:
				new_seq_name = new_seq_name + "ERRmidBC"
				# Do something here #
	
	# Save those without barcode
	seq_withoutbc_list = list(set(list(SeqDict.keys())) - set(seq_withbc_list))	
	for seq_withoutbc in seq_withoutbc_list:
		bc_leftover.write('>' + str(seq_withoutbc) + '\n' + str(SeqDict[seq_withoutbc].seq) + '\n')

	bc_blast.close()
	bc_trimmed.close()
	bc_leftover.close()
	bc_onlyF.close()
	bc_onlyR.close()
	bc_toomany.close()
	bc_invalid.close()

def DeBarcoder_SWalign(SeqDict, barcode_seq_filename, Output_folder, Output_prefix, search_range=25):
	"""Use Smith-Waterman local alignment to find barcodes, slow"""
	sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))

	barcode_seq = parse_fasta(barcode_seq_filename)
	out_stat = open(Output_folder + '/SWalign_barcode_out.txt', 'w')
	bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'a') # For writing the de-barcoded sequences

	count_matches = 0
	for each_rec in sorted(SeqDict):
		seq_to_search_F = str(SeqDict[each_rec].seq)[:search_range]
		seq_to_search_R = ReverseComplement(str(SeqDict[each_rec].seq))[:search_range]
		Match = None
		# Go through each barcode, do alignment, and see if there is a match or not
		for each_bc in barcode_seq:
			# Forward
			aln = sw.align(each_bc.seq, seq_to_search_F)
			aln_len = aln.r_end - aln.r_pos
			mismatches = len(each_bc.seq) - aln_len + aln.mismatches
			if mismatches <= 3:
				if not Match or aln.score > Match[1].score:
					Match = (each_bc.id, aln, mismatches, '+', each_bc.id)
			# Reverse
			aln = sw.align(each_bc.seq, seq_to_search_R)
			aln_len = aln.r_end - aln.r_pos
			mismatches = len(each_bc.seq) - aln_len + aln.mismatches
			if mismatches <= 3:
				if not Match or aln.score > Match[1].score:
					Match = (each_bc.id, aln, mismatches, '-', each_bc.id)

		if Match:
			#print Match[3], Match[1].dump()
			count_matches = count_matches + 1
			out_stat.write(str(each_rec) + '\t' + Match[0] + '\t' + str(Match[1].q_end) + '\t' + str(Match[2]) + '\t' + Match[3] + '\n')
			new_seq_name = Match[4] + '|' + str(each_rec) + '_recycled'
			if Match[3] == '+':
				bc_trimmed.write('>' + new_seq_name + '\n' + str(SeqDict[each_rec].seq)[Match[1].q_end:] + '\n')
			else:
				bc_trimmed.write('>' + new_seq_name + '\n' + ReverseComplement(str(SeqDict[each_rec].seq)[Match[1].q_end:]) + '\n')
	return count_matches

def dePrimer(FRprims, InFile, OutFile):
	'''Attempt to use Usearch to remove primers; not useful as it does not allow indels; stick to cutadapt'''
	#usearch -search_pcr greengenes.fa -db 16s_primers.fa -strand both -maxdiffs 3 -minamp 30 -maxamp 2000 -pcrout hits.txt -ampout amplicons.fasta
	usearch_cline = "%s -search_pcr %s -db %s -strand both -maxdiffs 3 -minamp 30 -maxamp 2000 -pcrout hits.txt -ampout %s" % (Usearch, InFile, FRprims, OutFile)

def doCutAdapt(Fprims, Rprims, InFile, OutFile, minimum_len=50): 
	'''This function removes the primers using the Cutadapt program. Replaces the DePrimer function'''
	
	# Build the forward and reverse primer portions of the cutadapt command lines by adding each primer in turn
	F_cutadapt_cline = ''
	for each_primer in Fprims: 
		each_primer = each_primer.upper() #in case the primers aren't uppercase (lower case nucs will confuse ReverseComplement)
		# F_cutadapt_cline = F_cutadapt_cline + '-g ' + each_primer + ' ' # F_cutadapt_cline as '-g GGACCTGGSCTYGCTGARGAGTG -g TCTGCMCATGCMATTGAAAGAGAG -g GAGYGTTTGGAATGTYTCWTTCCTYGG'
		F_cutadapt_cline = F_cutadapt_cline + '-g ^' + each_primer + ' ' 
		# e.g.,  F_cutadapt_cline as '-g ^GGACCTGGSCTYGCTGARGAGTG -g ^TCTGCMCATGCMATTGAAAGAGAG -g ^GAGYGTTTGGAATGTYTCWTTCCTYGG'
		""" The version with the "^" included anchors the fprimer so that it has to occur at the beginning of the sequence.
		This has the advantage of not picking up primers in the middle (in which case the first part--presumably including the 
		original barcode--would be erased). But it may leave behind too many primers (primers have to be within the error rate 
		of full-length) which may mess up the clustering. And R primers can still be found in the middle of the sequence. """
	
	R_cutadapt_cline = ''
	for each_primer in Rprims:
		each_primer = each_primer.upper() #in case the primers aren't uppercase (lower case nucs will confuse ReverseComplement)
		R_cutadapt_cline = R_cutadapt_cline + '-a ' + ReverseComplement(each_primer) + ' '

	# Build the complete cutadapt command line
	cutadapt_cline = '%s -O 15 -e 0.25 -n 2 -m %s %s %s %s > %s' % (Cutadapt, minimum_len, F_cutadapt_cline, R_cutadapt_cline, InFile, OutFile)
	# -O: Minimum overlap length. If the overlap between the read and the primer is shorter than -O, the read is not modified.
	# -e: Maximum allowed error rate (no. of errors divided by the length of the matching region)
	# -n: Try to remove primers at most -n times.
	
	process = subprocess.Popen(cutadapt_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate()
	if verbose_level in [1, 2]:
		print '**Primer-trimming results**\n'
		print err

def SplitBy(annotd_seqs_file, split_by = "locus-taxon", Multiplex_perBC_flag=True): #I can get rid of Output_folder? Output_folder,
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

def makeMapDict(mapping_file, locus, Multiplex_perBC_flag=True, DualBC_flag=False): #danger? Locus used somewhere else?
	map = open(mapping_file, 'rU') #open the correct mapping file
	outputfile_name = '_' + str(locus) + '.txt'

	#constructing a dict that maps each barcode/locus combination to a specific accession
	if Multiplex_perBC_flag and not DualBC_flag: # when there are multiple individuals sharing the same barcodes
		MapDict = {}
		for each_rec in map:
			each_rec = each_rec.strip('\n')
			map_barcode_name = each_rec.split('\t')[0]
			map_group_name = each_rec.split('\t')[1].upper()
			map_taxon_name = each_rec.split('\t')[2]
			key = map_barcode_name + '_' + map_group_name #BC01_A, BC01_B, BC010_C...
			MapDict[key] = map_taxon_name
	elif not Multiplex_perBC_flag and DualBC_flag:
		MapDict = {}
		for each_rec in map:
			each_rec = each_rec.strip('\n')
			map_barcodeF_name = each_rec.split('\t')[0]
			map_barcodeR_name = each_rec.split('\t')[1]
			map_taxon_name = each_rec.split('\t')[2]
			key = map_barcodeF_name + '^' + map_barcodeR_name #BCF01^BCR01, BCF02^BCR02, ...
			MapDict[key] = map_taxon_name
	elif not Multiplex_perBC_flag and not DualBC_flag: # when a barcode points to one individual; no multiplex per barcode; ignore "group" info and reads a different mapping file format
		MapDict = {}
		for each_rec in map:
			each_rec = each_rec.strip('\n')
			map_barcode_name = each_rec.split('\t')[0]
			map_taxon_name = each_rec.split('\t')[1]
			key = map_barcode_name #BC01_A, BC01_B, BC010_C...
			MapDict[key] = map_taxon_name
	
	map.close() #think this belongs in here now
	return MapDict

def annotateIt(filetoannotate, outFile, failsFile, Multiplex_perBC_flag=True, DualBC_flag=False, verbose_level=0):	
	"""Uses the blast results (against the reference sequence database) to assign locus and taxon, and write sequences 
	for a particular locus as specified by map_locus; returns a dictionary containing taxon-locus seq counts"""		
	
	# Blasts each sequence in the file (e.g., BC01.fa) against the reference sequences
	BlastSeq(filetoannotate, Output_folder + '/blast_refseq_out.txt', BLAST_DBs_folder + '/' + refseq_databasefile, evalue=0.0000001, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')
	
	SeqDict = SeqIO.index(filetoannotate, 'fasta') # Reads the sequences as a dict

	refseq_blast = open(Output_folder + '/blast_refseq_out.txt', 'rU') #read the ref_seq blast result
	groupsList = []
	locusList = []
	LocusTaxonCountDict = {}
	
	# Using blast matches to the reference sequences, and barcode <-> taxon mapping files, to assign 
	# each seq to a particular locus and taxon
	dictOfMapDicts = {} # A dictionary to store all of the map dictionaries
	for each_file, each_locus in zip(mapping_file_list, locus_list): 
		dictOfMapDicts[each_locus.upper()] = makeMapDict(each_file, each_locus, Multiplex_perBC_flag, DualBC_flag) # Note the Multiplex_perBC and DualBC flags (as True/False)

	annotated_seqs = open(outFile, "w")
	no_matches = open(failsFile, "w")

	#if verbose_level in [1,2]:
	#	sys.stderr.write("Annotating " + str(len(SeqDict)) + " records.\n")

	#count = 0
	#previous_seq_name = ''
	seq_processed_list = []
	for each_rec in refseq_blast:
		
		#count += 1
		#if verbose_level in [1,2]:
		#	if count == int(len(SeqDict)/4):
		#		sys.stderr.write( "One quarter done\n" ) # FWL; TODO; how to get this to print as the function proceeds, rather than just all at the end?
		#	elif count == int(len(SeqDict)/2):
		#		sys.stderr.write( "Half done\n" )
		#	elif count == int(len(SeqDict)*.75):
		#		sys.stderr.write( "Three quarters done\n" )

		each_rec = each_rec.strip('\n')	
		seq_name = each_rec.split('\t')[0] # The un-annotated sequence name, e.g., "BC02|m131213_174801_42153_c100618932550000001823119607181400_s1_p0/282/ccs;ee=7.2;"
		refseq_name = each_rec.split('\t')[1].replace(' ','') # The best-hit reference sequence name, e.g., "PGI_grC__C_diapA_BC17"
				
		# Get the key for retrieving taxon_name in dictOfMapDicts[locus_name]
		if Multiplex_perBC_flag:
			group_name = refseq_name.split('_')[1][-1:].upper() # Trying to grab the last letter of the second "word", i.e., the "A" in "grA"		
			locus_name = refseq_name.split('_')[0].upper() # The names are in the format ">ApP_grA__otherstuff". I.e., >Locus_group_otherstuff
			key = seq_name.split('|')[0] + '_' + group_name # Grabbing the barcode from the source seq, and the group from the matching ref seq.
			#i.e., gets the unique identifier that can link to a specific sample; i.e. BC01_A, BC01_B, BC01_C...
			if not group_name in groupsList: #keeping track of which groups are found, as a way of potentially diagnosing errors
				groupsList.append(group_name)
			if not locus_name in locusList: #keeping track of which loci are found, as a way of potentially diagnosing errors
				locusList.append(locus_name)	
		else:
			key = seq_name.split('|')[0] # Grabbing the barcode from the source seq
			locus_name = locus_list[0].upper()
			#i.e., gets the unique barcode that can link to a specific sample; i.e. BC01, BC02, BC03...		
		
		try: #use try/except to avoid the error when the key is not present in MapDict				
			taxon_name = dictOfMapDicts[locus_name][key] 
			#getting to the dict corresponding to this locus, and then finding that taxon that matches the barcode+group (the key)
			
			if Multiplex_perBC_flag:
				new_seq_name = taxon_name + '|' + locus_name + '|' + group_name + '|' + seq_name.replace(seq_name_toErase, '')
			else:
				new_seq_name = taxon_name + '|' + locus_name + '|' + seq_name.replace(seq_name_toErase, '')
				
			if seq_name not in seq_processed_list:
				annotated_seqs.write('>' + new_seq_name + '\n' + str(SeqDict[seq_name].seq) + '\n')
				try:
					LocusTaxonCountDict[taxon_name, locus_name] += 1 #as {('C_mem_6732', 'PGI'): 2, ('C_mem_6732', 'IBR'): 4} for example
				except:
					LocusTaxonCountDict[taxon_name, locus_name] = 1 #initiate the key and give count = 1
				seq_processed_list.append(seq_name)
		except:
			print "The combo", key, "wasn't found in", locus_name
			#print "(currently trying to find a match for", seq_name, ")\n"
			if Multiplex_perBC_flag:
				new_seq_name = locus_name + '|' + group_name + '|' + seq_name.replace(seq_name_toErase, '')
			else:
				new_seq_name = locus_name + '|' + seq_name.replace(seq_name_toErase, '')
			no_matches.write('>' + new_seq_name + '\n' + str(SeqDict[seq_name].seq) + '\n')
			seq_processed_list.append(seq_name)
			continue
	
	seq_no_hit = list(set(SeqDict.keys()) - set(seq_processed_list))
	print len(seq_no_hit)
	for each_rec in seq_no_hit:
		no_matches.write('>' + each_rec + '\n' + str(SeqDict[each_rec].seq) + '\n')
	
	refseq_blast.close()
	annotated_seqs.close()
	no_matches.close()

	if verbose_level in [1, 2]:
		print "The groups found are ", ', '.join(groupsList), "\nAnd the loci found are ", ', '.join(locusList), "\n"
	
	return LocusTaxonCountDict #as {('C_mem_6732', 'PGI'): 2, ('C_mem_6732', 'IBR'): 4} for example

def sortIt_length(file, verbose_level=0):
	"""Sorts clusters by seq length"""
	outFile = re.sub(r"(.*)\..*", r"\1_Sl.fa", file) # Make the outfile name by cutting off the extension of the infile name, and adding "_S1.fa"
	usearch_cline = "%s -sortbylength %s -fastaout %s" %(Usearch, file, outFile)
	#usearch_cline = "%s -sortbylength %s -output %s" %(Usearch, file, outFile) # Usearch 7 
	process = subprocess.Popen(usearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		print '\n**Usearch-sorting output on', file, '**\n'
		print err
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
		print '\n**Usearch-clustering output on', file, '**\n'
		print err
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
		print '\n**Uchime output on', file, '**\n'
		print err
	return outFile
	
'''TODO Add in a "chimera_mode" option that can be set to denovo or reference-based. 
Haven"t looked into how to implement this yet '''

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
		print '\n**Usearch-sorting output on', file, '**\n'
		print err
	return outFile

def muscleIt(file, verbose_level=0):
	"""Aligns the sequences using MUSCLE"""
	outFile = re.sub(r"(.*)\..*", r"\1.afa", file) # The rs indicate "raw" and thus python's escaping gets turned off
	muscle_cline = '%s -in %s -out %s' % (Muscle, file, outFile)
	process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout 
	if verbose_level == 2:
		print '\n**Muscle output on', file, '**\n'
		print err
	return outFile




################################################ Setup ################################################

if len(sys.argv) < 2:
	sys.exit(usage)
	
elif sys.argv[1] in ['-help', '-h', '-citation']:
	sys.exit(usage + citation)

else:
	print "PPP called with: \n\t", sys.argv, "\n" # Why does this print command (and the others, below)
	# not get printed earlier if using a redirect? ie, it says "Reading sequences.." but doesn't say "PPP called with:.."
	try:
		configuration = open(sys.argv[1], 'rU')
	except:
		sys.stderr.write('Error: Cannot open the configuration file\n')
		sys.exit(usage)

	## Setting up defaults ##
	mode = 0
	Check_chimeras = False
	Multiplex_per_barcode = False
	Dual_barcode = False
	Search_ends_only = False
	Recycle_bc = True
	Align = 0
	verbose_level = 0

	## Read-in the parameters and settings ##
	for line in configuration:
		line = line.strip(' \t\n\r').replace(' ', '').replace('\t', '').split('#')[0]
		setting_line = line.find('=')
		if setting_line != -1:
			setting = line.split('=')
			setting_name = setting[0]
			setting_argument = setting[1]
			if setting_name == 'Mode':
				mode = int(setting_argument)
			elif setting_name == 'Input_sequence_file':
				raw_sequences = setting_argument
			elif setting_name == "Align":
				Align = int(setting_argument)
			elif setting_name == 'Output_prefix':
				Output_prefix = setting_argument
			elif setting_name == 'Output_folder':
				Output_folder = setting_argument
				BLAST_DBs_folder = Output_folder + '_BlastDBs'
			elif setting_name == 'Barcode_blastDB':
				barcode_databasefile = setting_argument		
			elif setting_name == 'RefSeq_blastDB':
				refseq_databasefile = setting_argument
			elif setting_name == 'Locus_name':
				locus_list = setting_argument.upper().replace(' ', '').replace('\t', '').split(',') #needs the upper() now that LocusTaxonCountDict_unclustd has the loci in uppercase
			elif setting_name == 'Locus-barcode-taxon_map':
				mapping_file_list = setting_argument.replace(' ', '').replace('\t', '').split(',')
			elif setting_name == 'Usearch':
				Usearch = setting_argument
			elif setting_name == 'Cutadapt':
				Cutadapt = setting_argument
			elif setting_name == 'Muscle':
				Muscle = setting_argument
			elif setting_name == 'clustID':
				clustID = float(setting_argument)
			elif setting_name == 'clustID2':
				clustID2 = float(setting_argument)
			elif setting_name == 'clustID3':
				clustID3 = float(setting_argument)													
			elif setting_name == 'sizeThreshold':
				sizeThreshold = float(setting_argument)	
			elif setting_name == 'sizeThreshold2':
				sizeThreshold2 = float(setting_argument)
			elif setting_name == 'Forward_primer':
				Forward_primer = setting_argument.replace(' ', '').replace('\t', '').split(',')
			elif setting_name == 'Reverse_primer':
				Reverse_primer = setting_argument.replace(' ', '').replace('\t', '').split(',')					
			elif setting_name == 'seq_name_toErase':
				seq_name_toErase = setting_argument
			elif setting_name == 'verbose_level':
				verbose_level = int(setting_argument)	
			elif setting_name == 'split_type':
				split_type = setting_argument				
			elif setting_name == 'in_Barcode_seq_file':	
				barcode_seq_filename = setting_argument
			elif setting_name == 'in_RefSeq_seq_file':	
				refseq_filename = setting_argument
			elif setting_name == 'Dual_barcode':
				if setting_argument == '0':
					Dual_barcode = False
				elif setting_argument == '1':
					Dual_barcode = True
				else:
					sys.exit('Error: incorrect setting of Dual_barcode')			
			elif setting_name == 'Multiplex_per_barcode':	
				if setting_argument == '0':
					Multiplex_per_barcode = False
				elif setting_argument == '1':
					Multiplex_per_barcode = True
				else:
					sys.exit('Error: incorrect setting of Multiplex_per_barcode')	
			elif setting_name == 'Barcode_detection':
				if setting_argument == '0':
					Search_ends_only = False
				elif setting_argument == '1':
					Search_ends_only = True
				else:
					sys.exit('Error: incorrect setting of Barcode_detection')						
			elif setting_name == 'Recycle_no_barcoded_seq':
				if setting_argument == '0':
					Recycle_bc = False
				elif setting_argument == '1':
					Recycle_bc = True
				else:
					sys.exit('Error: incorrect setting of Recycle_no_barcoded_seq')
			elif setting_name == 'Recycle_chimeric_seq':
				if setting_argument == '0':
					Recycle_chimera = False
				elif setting_argument == '1':
					Recycle_chimera = True
				else:
					sys.exit('Error: incorrect setting of Recycle_chimeric_seq')

	if Recycle_bc:
		import swalign #so that users don't need to have this if not using this functionality

	print "Settings for this run:"
	print "Sequence file:\t", raw_sequences, "\n Loci:\t", locus_list

	if mode in [0,1]:
		print "Mapping files:\t", mapping_file_list
		print "Forward primers:\t", Forward_primer
		print "Reverse primers:\t", Reverse_primer
		if Dual_barcode:
			print "Expecting barcodes at each end of the sequence"
		else:
			print "Expecting barcodes on forward primers only"
		if Multiplex_per_barcode:
			print "Expecting barcodes to be shared across multiple taxa (genera, etc)"
		else:
			print "Expecting each barcode to be used for only a single taxon"
		if Search_ends_only:
			print "Primers will be looked for in the terminal XXXXX bases of each sequence, only"
		else:
			print "The full sequence will be searched for primers; internal primers may be found"
		if Recycle_bc:
			print '''If the BLAST-based approach doesn't find a barcode in a particular sequence
			the sequence will be re-searched using a Smith-Waterman pairwise alignment approach'''
		else: 
			print "Primer-searching will be done using the BLAST-based approach, only"
		if Recycle_chimera:
			print "Inter-locus chimeras will be split into their component pieces and fed into the pipeline"
		else:
			print "Inter-locus chimeras will be removed from the analysis"
	if mode != 2: # ie, some clustering will be done
		print "Similarity cut-off for clustering:\t", clustID, "\t", clustID2, "\t", clustID3
		print "Cluster size for retention:\t", sizeThreshold, "\t", sizeThreshold, "\t", sizeThreshold2
		# may want to add more outputs here
	




################################################ RUN YO!!! ########################################

if mode in [0,1]: # Make blast databases and read the raw sequences
	if os.path.exists(BLAST_DBs_folder): # overwrite existing folder
		shutil.rmtree(BLAST_DBs_folder)
	os.makedirs(BLAST_DBs_folder)
	os.chdir(BLAST_DBs_folder)
	makeBlastDB('../' + refseq_filename, refseq_databasefile) # and one of the reference sequences
	makeBlastDB('../' + barcode_seq_filename, barcode_databasefile) # one of the barcodes
	#if Dual_barcode:
		#makeBlastDB('../' + barcode_seq_filename2, barcode_databasefile + '2')
	os.chdir('..')
	
	## Read sequences ##
	sys.stderr.write('Reading sequences...\n')
	SeqDict = SeqIO.index(raw_sequences, 'fasta') # Read in the raw sequences as dictionary, using biopython's function
	count_total_input_sequences = len(SeqDict)
	sys.stderr.write('\t' + str(count_total_input_sequences) + ' sequences read\n')

	## Make output folder ##
	if os.path.exists(Output_folder): # overwrite existing folder
		shutil.rmtree(Output_folder)
	os.makedirs(Output_folder)

if mode == 0: # QC mode
	## Check chimeras ##
	sys.stderr.write('Checking for inter-locus chimeric sequences (concatemers)...\n')
	if not Recycle_chimera:
		chimeras_file = Output_prefix + '_0_chimeras.fa'
		non_chimeras_file = Output_prefix + '_0_nonchimeras.fa'
		chimera_dict = CheckChimericLoci(raw_sequences, Output_folder + '/' + 'blast_full_refseq_out.txt', Output_folder + '/' + non_chimeras_file, Output_folder + '/' + chimeras_file, BLAST_DBs_folder + '/' + refseq_databasefile, SeqDict, SplitChimera=False)
	else:
		chimeras_file = Output_prefix + '_0_chimeras.fa'
		non_chimeras_file = Output_prefix + '_0_nonchimeras+split.fa'
		chimera_dict = CheckChimericLoci(raw_sequences, Output_folder + '/' + 'blast_full_refseq_out.txt', Output_folder + '/' + non_chimeras_file, Output_folder + '/' + chimeras_file, BLAST_DBs_folder + '/' + refseq_databasefile, SeqDict, SplitChimera=True)
	
	count_chimeric_sequences = len(chimera_dict)
	raw_sequences = Output_folder + '/' + non_chimeras_file
	SeqDict = SeqIO.index(raw_sequences, 'fasta')
	sys.stderr.write('\t' + str(count_chimeric_sequences) + ' inter-locus chimeric sequences found\n')

if mode in [0,1]: # Run the full annotating, clustering, etc.	
	## Remove barcodes ##
	if Dual_barcode:
		sys.stderr.write('Removing dual barcodes...\n')
		DeBarcoder_dual(raw_sequences, BLAST_DBs_folder + '/' + barcode_databasefile, SeqDict)
	else:
		sys.stderr.write('Removing barcodes...\n')
		if not Search_ends_only:
			DeBarcoder(raw_sequences, BLAST_DBs_folder + '/' + barcode_databasefile, SeqDict, Output_folder, Output_prefix)
		else:
			DeBarcoder_ends(SeqDict, BLAST_DBs_folder + '/' + barcode_databasefile, Output_folder, Output_prefix, search_range=25)
	count_seq_w_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa')
	count_seq_wo_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa')
	count_seq_w_toomany_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa')
	sys.stderr.write('\t' + str(count_seq_w_bc) + ' sequences have barcode\n')
	sys.stderr.write('\t' + str(count_seq_wo_bc) + ' sequences have no barcode\n')
	sys.stderr.write('\t' + str(count_seq_w_toomany_bc) + ' sequences have too many barcodes\n')

	if Recycle_bc:
		sys.stderr.write('Looking for primers in the no-primer sequences, using Smith-Waterman pairwise alignment...\n')
		SeqDict_no_bc = SeqIO.index(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'fasta') # Read in the raw sequences as dictionary, using biopython's function
		count_seq_recycled = DeBarcoder_SWalign(SeqDict_no_bc, barcode_seq_filename, Output_folder, Output_prefix, search_range=25)
		sys.stderr.write('\t' + str(count_seq_recycled) + ' sequences recycled from ' + str(count_seq_wo_bc) + ' sequences\n')
	
	## Remove primers ##
	sys.stderr.write('Removing primers...\n')
	primer_trimmed_file = Output_prefix + '_2_pr_trimmed.fa'
	doCutAdapt(Fprims = Forward_primer, Rprims = Reverse_primer, InFile = Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', OutFile = Output_folder + '/' + primer_trimmed_file)
	count_seq_pr_trimmed = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_2_pr_trimmed.fa')
	sys.stderr.write('\t' + str(count_seq_pr_trimmed) + ' sequences survived after primer-trimming\n')

	## Annotate the sequences with the taxon and locus names, based on the reference sequences ##
	sys.stderr.write('Annotating seqs...\n')
	toAnnotate = primer_trimmed_file
	annoFileName = Output_prefix + '_3_annotated.fa'
	LocusTaxonCountDict_unclustd = annotateIt(filetoannotate = Output_folder + '/' + toAnnotate, outFile = Output_folder + '/' + annoFileName, Multiplex_perBC_flag = Multiplex_per_barcode, DualBC_flag = Dual_barcode, failsFile = Output_folder + '/' + Output_prefix + '_3_unclassifiable.fa', verbose_level = verbose_level)
	count_seq_annotated = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_3_annotated.fa')
	count_seq_unclassifiable = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_3_unclassifiable.fa')
	sys.stderr.write('\t' + str(count_seq_annotated) + ' sequences annotated\n')
	sys.stderr.write('\t' + str(count_seq_unclassifiable) + ' sequences cannot be classified\n')

	## Move into the designated output folder ##
	os.chdir(Output_folder)

	## Split sequences into separate files/folders for each locus ##
	sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
	locusCounts = SplitBy(annotd_seqs_file = annoFileName, split_by = "locus", Multiplex_perBC_flag = Multiplex_per_barcode) 

	## Split the locus files by barcode, and cluster each of the resulting single locus/barcode files
	sys.stderr.write('Clustering/dechimera-izing seqs...\n')
	all_folders_loci = locusCounts.keys() # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
		# folders) and they correspond to the counts for each.
	LocusTaxonCountDict_clustd = {}

	for each_folder in all_folders_loci: 
		os.chdir(each_folder)
		sys.stderr.write('\nWorking on: ' + each_folder + '...\n')
		print '\nWorking on ', each_folder, ' ...\n'
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
				print "Working on ", bcode_folder
				os.chdir(bcode_folder)
				
				print "\tFirst clustering"
				print "\nAttempting to sort: ", bcode_folder + ".fa"

				sorted_length = sortIt_length(file = bcode_folder + ".fa", verbose_level = verbose_level)
				clustered1 = clusterIt(file = sorted_length, clustID = clustID, round = 1, verbose_level = verbose_level)

				print "\tFirst chimera slaying expedition"
				deChimered1 = deChimeIt(file = clustered1, round = 1, verbose_level = verbose_level)
				
				print "\tSecond clustering"
				sorted_size1 = sortIt_size(file = deChimered1, thresh = sizeThreshold, round = 1, verbose_level = verbose_level)
				clustered2 = clusterIt(file = sorted_size1, clustID = clustID2, round = 2, verbose_level = verbose_level)
				
				print "\tSecond chimera slaying expedition"
				deChimered2 = deChimeIt(file = clustered2, round = 2, verbose_level = verbose_level)
				
				print "\tThird clustering"
				sorted_size2 = sortIt_size(file = deChimered2, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
				clustered3 = clusterIt(file = sorted_size2, clustID = clustID3, round = 3, verbose_level = verbose_level)
				
				print "\tThird chimera slaying expedition\n"
				deChimered3 = deChimeIt(file = clustered3, round = 3, verbose_level = verbose_level)

				# Why are we sorting again? I guess this gives us the chance to remove clusters smaller than sizeThreshold2
				sorted_size3 = sortIt_size(file = deChimered3, thresh = sizeThreshold2, round = 3, verbose_level = verbose_level)

				clustered_seq_file = parse_fasta(sorted_size3)
				for each_seq in clustered_seq_file:
					taxon_name = str(each_seq.id).split('|')[0].split('=')[-1] # for example, get C_dia_5316 from centroid=centroid=C_dia_5316|ApP|C|BC02|_p0/158510/ccs;ee=1.9;;seqs=6;seqs=18;size=27;
					
					try:
						LocusTaxonCountDict_clustd[taxon_name, each_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
						# The locus names are the same as each_folder
					except:
						LocusTaxonCountDict_clustd[taxon_name, each_folder] = 1		
				os.chdir("..") # To get out of the current barcode folder and ready for the next one
		os.chdir("..") # To get out of the current locus folder and ready for the next one


	## Put all the sequences together ##
	sys.stderr.write('Putting all the sequences together......\n\n')
	for each_folder in all_folders_loci: # Looping through each of the locus folders
		outputfile_name = '_' + str(each_folder) + '_clustered.txt' # "each_folder" is also the name of the locus
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
		
	## Clean-up the sequence names ##
	print "Cleaning up the file names"
	sys.stderr.write("Cleaning up the file names\n")

	fastas = glob.glob("_*.txt")
	for file in fastas:
		print "Cleaning up the sequence names in " + file + "\n"
		fasta_cleaned = open(str(file).replace(".txt", "_renamed.fa"), 'w') # Make a new file with .fa instead of .txt
		parsed = parse_fasta(file)
		for seq in parsed:
			seq.id = re.sub(r"seqs=\d*", r"", seq.id)
			seq.id = re.sub(r"ccs.ee=[\d\.]*", r"", seq.id)
			seq.id = seq.id.replace("centroid=", "").replace(";", "").replace("/","_")
			fasta_cleaned.write(str('>' + seq.id + "\n" + seq.seq + "\n"))
		fasta_cleaned.close()


if mode == 2: # Just split the seqs
	annoFileName = Output_prefix + '_3_annotated.fa'
	print "The file to be split is ", annoFileName
	print "And it should be split by ", split_type
	print "And it should be in ", Output_folder

	#os.mkdir(Output_folder)
	os.chdir(Output_folder)
	splits_count_dict = SplitBy(annotd_seqs_file = annoFileName, split_by = split_type, Multiplex_perBC_flag = Multiplex_per_barcode)
	os.chdir("..")
	if verbose_level == 3:
		print "Counts: "
		for each_split in splits_count_dict:
			print each_split, '\t', splits_count_dict[each_split]

if mode == 3: # Just do some clustering
	os.chdir(Output_folder)
	tocluster = glob.glob("*") # Will attempt to cluster all files (and folders!) in the Output_folder directory
	print "going to try to cluster these files: ", tocluster
	for thisfile in tocluster:
				
		print "\nAttempting to sort sequences in: ", thisfile
		print "\tFirst clustering"

		outFile = sortIt_length(file = thisfile, verbose_level = verbose_level)
		outFile = clusterIt(file = outFile, clustID = clustID, round = 1, verbose_level = 0)
		
		# print "Before first chimera hunt, outFile is", str(outFile)

		print "\tFirst chimera slaying expedition"
		outFile = deChimeIt(file = outFile, round = 1, verbose_level = 0)
		
		print "\tSecond clustering"
		outFile = sortIt_size(file = outFile, thresh = sizeThreshold, round = 1, verbose_level = 0)
		outFile = clusterIt(file = outFile, clustID = clustID2, round = 2, verbose_level = 0)
		
		print "\tSecond chimera slaying expedition"
		outFile = deChimeIt(file = outFile, round = 2, verbose_level = 0) 
		
		print "\tThird clustering"
		outFile = sortIt_size(file = outFile, thresh = sizeThreshold, round = 2, verbose_level = 0)
		outFile = clusterIt(file = outFile, clustID = clustID3, round = 3, verbose_level = 0)
		
		print "\tThird chimera slaying expedition\n"
		outFile = deChimeIt(file = outFile, round = 3, verbose_level = 0) 

		outFile = sortIt_size(file = outFile, thresh = sizeThreshold2, round = 3, verbose_level = 0)

	## Cleaning seqs and getting them into locus-specific files ##
	tosummarize = glob.glob("*Ss3.fa")
	LocusTaxonCountDict_clustd = {}
	
	loci = []
	for each_file in tosummarize:
		thisfasta = parse_fasta(each_file)
		# Making the output file name from the first four underscore-separated elements of the original file name, 
		# and sticking them back together with underscores
		outFileName = "_".join(each_file.split("_")[0:4]) + "_clean.fa"
		fasta_cleaned = open(outFileName, 'w') 
		for each_seq in thisfasta:
			taxon_name = str(each_seq.id).split('|')[0].split('=')[-1] # for example, get C_dia_5316 from centroid=centroid=C_dia_5316|ApP|C|BC02|_p0/158510/ccs;ee=1.9;;seqs=6;seqs=18;size=27;	
			locus_name = str(each_seq.id).split("|")[1]
			
			if not locus_name in loci: # storing the loci that are present in a list, so that I can loop over it later
				loci.append(locus_name)

			try:
				LocusTaxonCountDict_clustd[taxon_name, locus_name] += 1  # {('C_dia_5316', 'ApP'): 28} for example
			except:
				LocusTaxonCountDict_clustd[taxon_name, locus_name] = 1

			each_seq.id = re.sub(r"seqs=\d*", r"", each_seq.id)
			each_seq.id = re.sub(r"ccs;ee=[\d\.]*", r"", each_seq.id)
			each_seq.id = each_seq.id.replace("centroid=", "").replace(";", "")
			fasta_cleaned.write(str('>' + each_seq.id + "\n" + each_seq.seq + "\n"))
		fasta_cleaned.close()

	## Put all the sequences together ##
	sys.stderr.write('Putting all the sequences together......\n\n')
	cleanfiles = glob.glob("*_clean.fa")
	outputfiles = {}

	for locus in loci: # created the "loci" list in the summarizing loop, above
		outputfile_name = "_" + str(locus) + ".txt"
		outputfiles[locus] = open(outputfile_name, 'w')
	for cleanfile in cleanfiles:
		#I thnk I've got this fixed now, but keeping these comments here for now in case I missed something
		# locus = str(cleanfile).split("_")[-2] # find out what locus the file contains
		
		# # This is messing up because the *_clean.fa file name, when clustering data that has been split by 
		# # group or barcode etc, doesn't contain the locus name
		# # TODO - fix this. What happens if I attempt to cluster a group that contains multiple loci?
		# shutil.copyfileobj(open(cleanfile,'rb'), outputfiles[locus]) # Add each file to the right output file

		sequences = parse_fasta(cleanfile)
		for eachSeq in sequences:
			locus = eachSeq.id.split("|")[1]
			SeqIO.write(eachSeq, outputfiles[locus], "fasta")

	for thisfile in outputfiles:
		outputfiles[thisfile].close()

	

#### Producing a summary #### 
if mode != 2: # Need to decide what sort of summary would be appropriate if just splitting
	count_output = open('counts.xls', 'w')
	count_output.write('Total input sequences:\t' + str(count_total_input_sequences) + '\n')
	if mode == 0:
		count_output.write('Concatemers (multi-locus seqs):\t' + str(count_chimeric_sequences) + '\n')
	if mode in [0,1]:
		count_output.write('Sequences with barcodes:\t' + str(count_seq_w_bc) + '\n')
		count_output.write('Sequences without barcodes:\t' + str(count_seq_wo_bc) + '\n')
		count_output.write('Sequences with too many barcodes:\t' + str(count_seq_w_toomany_bc) + '\n')
		count_output.write('Sequences annotated:\t' + str(count_seq_annotated) + '\n')
		count_output.write('Sequences that cannot be classified:\t' + str(count_seq_unclassifiable) + '\n')

	print "\n**Raw reads per accession per locus**"
	count_output.write('**Raw reads per accession per locus**\n')
	taxon_list =[]

# getting a list of all the taxa with sequences, from the count dictionary
# Using the counts of unclustered sequences so as to not miss any taxa
if mode in [0,1]:
	for i in range(0,len(LocusTaxonCountDict_unclustd)): 
	# The keys for this dictionary are two-part lists, e.g., [('C_mem_6732', 'IBR'), ('C_mem_6732', 'PGI'), ('C_dou_111', 'IBR')]
		if not LocusTaxonCountDict_unclustd.keys()[i][0] in taxon_list:
			taxon_list.append(LocusTaxonCountDict_unclustd.keys()[i][0])
if mode == 3:
	print
	# do something
	#TODO LocusTaxonCountDict_unclustd doesn't exist in mode 3, but we'd still like to summarize. How to do so?

print '\t', '\t'.join(locus_list)
count_output.write('\t' + '\t'.join(locus_list) + '\n')
for each_taxon in set(taxon_list): 
	print each_taxon, '\t',
	count_output.write(each_taxon + '\t')
	for each_locus in locus_list:
		try:
			print LocusTaxonCountDict_unclustd[each_taxon, each_locus], '\t', 
			count_output.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
		except:
			print '0', '\t', 
			count_output.write('0\t')
	print 
	count_output.write('\n')
print 
print "\n**Final clustered sequences per accession per locus**"
count_output.write('\n**Final clustered sequences per accession per locus**\n')
print '\t', '\t'.join(locus_list)
count_output.write('\t' + '\t'.join(locus_list) + '\n')	
for each_taxon in set(taxon_list):
	print each_taxon, '\t',
	count_output.write(each_taxon + '\t')		
	for each_locus in locus_list:
		try:
			print LocusTaxonCountDict_clustd[each_taxon, each_locus], '\t',
			count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
		except:
			print '0', '\t',
			count_output.write('0\t') 
	print
	count_output.write('\n')		
print 

print '\n**Allele/copy/cluster/whatever count by locus**'
count_output.write('\n**Allele/copy/cluster/whatever count by locus**\n')	
# I think this was breaking if a locus had no sequences, and thus that file is not created. Going to try "try"
for each_locus in locus_list:
	file_name = '_' + str(each_locus) + '_clustered.txt'
	try: #I'm hoping that this will stop the program from crashing if a locus has no sequences
		seq_no = len(parse_fasta(file_name))
		print '\t', each_locus, ':', seq_no
		count_output.write(str(each_locus) + '\t' + str(seq_no) + '\n')
	except:
		print '\t', each_locus, ':', 0
		count_output.write(str(each_locus) + '\t0\n')			
print

## Aligning the sequences ##
if Align: # Aligning can be turned on/off in the configuration file
	fastas = glob.glob("_*.fa")
	""" If you want to align the versions that haven't had their names cleaned up (ie, still have info
		on the expected number of errors, etc., change the glob to operate on _*.txt"""
	for file in fastas:
		print "Aligning ", file, "\n"
		sys.stderr.write("Aligning " + file + "\n")
		outFile = muscleIt(file, verbose_level)
	
# 	## Aligning the sequences ##
# if Align: # Aligning can be turned on/off in the configuration file
# 	fastas = glob.glob("_*")
# 	for file in fastas:
# 		print "Aligning ", file, "\n"
# 		sys.stderr.write("Aligning " + file + "\n")
# 		muscleIt(file, verbose_level)

	

