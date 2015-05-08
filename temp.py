#!/usr/bin/env python
import swalign
from Bio import SeqIO

sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))

SeqDict = SeqIO.index('PBR3_reduced_forTesting_seqname.fa', 'fasta') # Read in the raw sequences as dictionary, using biopython's function

def parse_fasta(infile):
	"""Reads in a fasta, returns a list of biopython seq objects"""
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq] 

def ReverseComplement(seq):
	"""Returns reverse complement sequence, ignores gaps"""
	seq = seq.replace(' ','') # Remove spaces
	seq = seq[::-1] # Reverse the sequence
	basecomplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y':'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W', 'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B'} # Make a dictionary for complement
	letters = list(seq) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters) 

def DeBarcoder_SWalign(SeqDict, barcode_seq_filename, outputfile_bc_stat, outputfile_bc_trimmed, search_range=25):
	sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))

	barcode_seq = parse_fasta(barcode_seq_filename)
	out_stat = open(outputfile_bc_stat, 'w')
	out_trimmed = open(outputfile_bc_trimmed, 'w')

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
			print Match[3], Match[1].dump()
			out_stat.write(str(each_rec) + '\t' + Match[0] + '\t' + str(Match[1].q_end) + '\t' + str(Match[2]) + '\t' + Match[3] + '\n')
			new_seq_name = Match[4] + '|' + str(each_rec)
			if Match[3] == '+':
				out_trimmed.write('>' + new_seq_name + '\n' + str(SeqDict[each_rec].seq)[Match[1].q_end:] + '\n')
			else:
				out_trimmed.write('>' + new_seq_name + '\n' + ReverseComplement(str(SeqDict[each_rec].seq)[Match[1].q_end:]) + '\n')
	return

DeBarcoder_SWalign(SeqDict, 'barcode_seq_plate1_2_Dep.txt', 'out', 'out.fas')

def _tag_aln_check(aln, seqlen, barcodelen, orientation, edit=0, pos=0):
	r_len = aln.r_end - aln.r_pos
	mismatches = barcodelen - r_len + aln.mismatches

	if mismatches <= edit:
		if orientation == '5' and aln.q_pos <= pos:  # allow at most {pos} extra bases at start
			return True, ''
		elif orientation == '3' and aln.q_end >= seqlen - pos:  # allow at most {pos} extra bases at end
			return True, ''
		else:
			return False, 'Not at right location'
	else:
		return False, 'too many mismatches'

	return False, '?'

def check_tags(barcodes, seq, edit, pos, allow_revcomp=False, verbose=False):
	'''
	For each barcode, pull out the appropriate 5' or 3' sub sequence from {seq}. Then
	run a local alignment of the barcode to the subseq. If a good match is found, return
	it; otherwise, find the best match and return that.
	returns a multi-tuple:
		valid, (tag, valid_seq, edits, reason_if_fail)
	For the alignments, the reference is the barcode, the query is the subset of the read
	that is possibly the barcode (5'/3' subseq)
	'''

	best = None


	for tag in barcodes:
		barcodeseq, orientation, strip = barcodes[tag]
		if orientation == '5':
			testseq = seq[:len(barcodeseq) + edit + pos]
		else:
			testseq = seq[-1 * (len(barcodeseq) + edit + pos):]

		aln = sw.align(barcodeseq, testseq)
		valid, reason = _tag_aln_check(aln, len(testseq), len(barcodeseq), orientation, edit, pos)
		if verbose:
			print 'Testing tag: %s vs %s' % (str(barcodes[tag]), testseq)
			aln.dump()
			print valid, reason
		if valid:
			return True, (tag, aln, True, '')

		if not best or aln.score > best[1].score:
			best = (tag, aln, True, reason)

		if allow_revcomp:
			if orientation == '5':
				testseq = seq[-1 * (len(barcodeseq) + edit + pos):]
			else:
				testseq = seq[:len(barcodeseq) + edit + pos]

			aln = sw.align(revcomp(barcodeseq), testseq)
			valid, reason = _tag_aln_check(aln, len(testseq), len(barcodeseq), '5' if orientation == '3' else '3', edit, pos)
			if verbose:
				print 'Testing tag: %s [rc] vs %s' % (str(barcodes[tag]), testseq)
				aln.dump()
				print valid, reason
			if valid:
				return True, (tag, aln, False, '')

			if not best or aln.score > best[1].score:
				best = (tag, aln, False, reason)

	if verbose:
		print 'BEST: ', best
		best[1].dump()

	return False, best