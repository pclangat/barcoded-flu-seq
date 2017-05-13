#!/usr/bin/python
## ~/scripts/primerID_pipeline.py <sample_name>

## Takes as input the sample name
## STEP ONE: 

'''Requirements: 
requires pip: brew install pip
requires biopython: pip install biopython
requires regex: pip install regex'''

##Import modules
import sys, os, regex, subprocess, operator
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio.Align import AlignInfo

###METHODS
def get_sample_reads(prfx):
	qc_paired_fastq = prfx+'.qc.fq'
	qc_paired_fasta = prfx+'.qc.fas'
	
	SeqIO.convert(qc_paired_fastq, "fastq", qc_paired_fasta, "fasta")
	reads = SeqIO.index(qc_paired_fasta, 'fasta')
	return reads

def check_barcodes(reads_dict):
	matches = []
	intact_barcoded_seqs = 0
	barcoded_seqs = defaultdict(list)
	
	##For each read, look for reverse primer (either full or partial match)
	for read in reads_dict:
		sequence = reads_dict[read].seq
		primer_pos, oriented_sequence, match_depth = check_rev_primer_match(rev_primer, sequence)
		matches.append(match_depth)
		
		##If rev primer sequence is found, check if intact rev barcode pattern exists
		if primer_pos > 0:
			rev_primer_seq = oriented_sequence[primer_pos:]
			intact_rev_barcode = regex.findall(pattern, str(rev_primer_seq))
			
			##If intact rev barcode pattern, add barcode and associated sequences to barcodes dict
			if intact_rev_barcode:
				intact_barcoded_seqs += 1
				barcode = (intact_rev_barcode[0]).split(rev_primer)[1]
				barcoded_seqs[barcode].append(oriented_sequence)
	return matches, intact_barcoded_seqs, barcoded_seqs
	
def check_rev_primer_match(r_primer, seq):
	match = 'none'
	pos = seq.find(r_primer)
	
	##if match found
	if pos > -1:
		match = 'full'
	##if no match found
	else:
		##try reverse complement of sequence
		rc_seq = seq.reverse_complement()
		pos = rc_seq.find(r_primer)
		
		##if match, re-orient sequence
		if pos > -1:
			seq = rc_seq
			match = 'full'
		else:
			##if still no match, try searching for partial match
			pos, seq, match = check_partial_rev_primer_match(r_primer, seq)
	return pos, seq, match
	
def check_partial_rev_primer_match(r_primer, seq):
	pos = -1
	match = 'none'
	
	pattern = '('+r_primer+'){e<=2}'
	fuzzy_match = regex.search(pattern, str(seq), regex.BESTMATCH)
	
	##if partial match found
	if fuzzy_match:
		match = 'partial'
		pos = fuzzy_match.span()[1]
	##if no partial match found
	else:
		##try reverse complement of sequence
		rc_seq = seq.reverse_complement()
		fuzz_match = regex.search(pattern, str(rc_seq), regex.BESTMATCH)
		
		##if partial match, re-orient sequence:
		if fuzzy_match:
			seq = rc_seq
			match = 'partial'
			pos = fuzzy_match.span()[1]
	return pos, seq, match
		
def fuzzy_substring(n, h):
	"""Calculates the fuzzy match of needle in haystack,
	using a modified version of the Levenshtein distance
	algorithm.
	The function is modified from the levenshtein function
	in the bktree module by Adam Hupp"""
	if len(n) > len(h):
		needle = h
		haystack = n
	else:
		needle = n
		haystack = h
		
	m, n = len(needle), len(haystack)
    # base cases
	if m == 1:
		return not needle in haystack
	if not n:
		return m
	row1 = [0] * (n+1)
	for i in range(0,m):
		row2 = [i+1]
		for j in range(0,n):
			cost = ( needle[i] != haystack[j] )
			row2.append( min(row1[j+1]+1, # deletion
                               row2[j]+1, #insertion
                               row1[j]+cost) #substitution
                           )
		row1 = row2
		#print row1
	#print min(row1)
	#print (0.3*m)
	#if min(row1) < (0.3*m):
	#	primer_pos = row1.index(min(row1)) - m
	#else:
	#	primer_pos = -1
	#return primer_pos
	return min(row1)
	
##Generate alignment for all seqs with same barcode
def get_alignment(records):
	##open up muscle
	child = subprocess.Popen(str(muscle_cline),
							stdin=subprocess.PIPE,
							stdout=subprocess.PIPE,
							stderr=subprocess.PIPE,
							universal_newlines=True,
							shell=(sys.platform!="win32"))
	SeqIO.write(records, child.stdin, "fasta")
	child.stdin.close()
	align = AlignIO.read(child.stdout, "fasta")
	return align

def generate_consensus(alignment):
	align_summary = AlignInfo.SummaryInfo(alignment)
	consensus = align_summary.dumb_consensus(threshold=0.51, ambiguous='N')
	return consensus
				
###INITIALISATIONS
##Paths to software
muscle_path = '/Users/pclangat/Software/muscle/muscle3.8.31_i86darwin32'

##These will become user-accepted later
min_length = 125
min_phred = 20
min_barcode_count = 2

fwd_primer = "CGGGGAAAATATGCAACAATCCT"
rev_primer = "GAGGGTTTCACTTGGACTGGG"
full_rev_primer = "GAGGGTTTCACTTGGACTGGGNNNNANNNNANNNNAAAGCAGTGGTATCAACGCA"
pattern = rev_primer+'[A-Z]{4}A[A-Z]{4}A[A-Z]{4}A'
#pattern = regex.compile(pattern)
pattern = regex.compile('('+pattern+'){e<=5}')

muscle_cline = MuscleCommandline(muscle_path)

reference_file = '/Users/pclangat/Desktop/Projects/2-PrimerID_pipeline/1-daniel_MiSeq_2017-03-01/reference.txt'
reference_seq = ''
with open(reference_file, 'r') as fh:
	for line in fh:
		reference_seq += line.strip()
#print reference_seq
		
###MAIN
if __name__ == '__main__':
	##STEP 1: Parse user input
	try:
		prefix = sys.argv[1]
	except:
		print('[USAGE]: %s sample_prefix' % sys.argv[0] )
		sys.exit(1)

	###STEP 2: Convert qc fastq files to fasta file and load records as dictionary index (i.e. does not save all into memory, good for large fastq files)
	reads_dict = get_sample_reads(prefix)
	print "\n>>>BARCODE FILTERING & MOLECULAR COUNTING SUMMARY<<<"		
	print "[INFO 1]: Total input sequences (after pairing and QC): %s" % len(reads_dict)
	
	###STEP 3: Find reverse primer (gene portion) in read sequence, count whether full or partial or no match
	matches, intact_barcoded_seqs, barcoded_seqs = check_barcodes(reads_dict)
	reads_dict.close()
	m = len(matches) - matches.count('none')
	percent_m = (m*100.0/len(matches))
	percent_i = (intact_barcoded_seqs*100.0/m)
	percent_b = (len(barcoded_seqs)*100.0/intact_barcoded_seqs)
	
	print "[INFO 2]: Sequences with reverse primer region: %s (%.1f%%)" % (m, percent_m)
	print "... Full matches to reverse primer: %s" % matches.count('full')
	print "... Partial matches to reverse primer: %s" % matches.count('partial')
	print "... Reverse primer not found: %s" % matches.count('none')
	print "[INFO 3]: Sequences with intact barcodes: %s (%.1f%%)" % (intact_barcoded_seqs, percent_i)
	print "... Representing # of unique barcodes: %s (%.1f%%)" % (len(barcoded_seqs), percent_b)

	###STEP 4: Count number of seqs associated to each barcode (i.e. barcode multiplicity/count), decide whether to process
	print "[INFO 4]: Checking barcode multiplicity and filtering"
	barcodes_to_process_count = 0
	barcode_distribution = {}
	consensus_records = []
	
	##For each barcode
	b=1
	for barcode in barcoded_seqs:
		print "barcode %s" % b
		barcode_count = len(barcoded_seqs[barcode])	
		
		##Add frequency of barcode counts to a dictionary for count distribution
		if barcode_count in barcode_distribution:
			barcode_distribution[barcode_count] += 1
		else:
			barcode_distribution[barcode_count] = 1
			
		##If barcode multiplicity is greater than minimum (simple case: 2)
		##add to list of records of barcoded seqs to process
		if barcode_count > min_barcode_count:
			barcodes_to_process_count += 1
			i = 1
			barcode_group_records = []
			
			##for each sequence in that barcode group, name it and create seq record
			for sequence in barcoded_seqs[barcode]:
				id_name = '%s-%s-%s' % (barcode, barcode_count, i)
				i += 1
				record = SeqRecord(sequence, id=id_name, description='')
				barcode_group_records.append(record)
			
			##Take barcode_group and align it
			print "getting alignment"
			aligned_barcode_group = get_alignment(barcode_group_records)
			
			##Get consensus of alignment
			consensus_name = '%s-%s-%s' % (prefix, barcode, barcode_count)
			consensus_seq = generate_consensus(aligned_barcode_group)

			if 'N' not in consensus_seq:
				consensus_rec = SeqRecord(consensus_seq, id=consensus_name, description='')
				consensus_records.append(consensus_rec)	
		b+=1
		
	print "... Barcode multiplicity (read count) distribution:"
	for i in barcode_distribution:
		print "\t%s\t%s" % (i, barcode_distribution[i])
	
	percent_p = 100.0*barcodes_to_process_count/len(barcoded_seqs)
	percent_c = 100.0*len(consensus_records)/barcodes_to_process_count
	
	print "... Barcodes with multiplicity >%s to process: %s (%.1f%%)" % (min_barcode_count, barcodes_to_process_count, percent_p)
	print "... Total unambiguous consensus (i.e. template) sequences: %s (%.1f%%)" % (len(consensus_records), percent_c)
	print "... Output generated to file %s.barcode_filtered.fas with headers '>Sample_Name-Barcode-#Reads'" % prefix
	##Put barcodes to process into a fasta file
	SeqIO.write(consensus_records, '%s.barcode_filtered.fas' % prefix, 'fasta')
	
	###STEP 5: Compare template sequences and group matching ones together and count #templates
	print "[INFO 5]: Checking template sequences for unique sequence variants"
	template_sequences = []
	unique_sequences = {}
	
	##Creates ordered list of pre-barcode sequences, sorted by sequence length
	for record in consensus_records:
		seq = record.seq
		seq = str(seq[:-34])
		template_sequences.append(seq)
	template_sequences.sort(key = len, reverse=True)
	
	##Check each template sequence for unique sequences
	for seq in template_sequences:
		is_substring = False
		
		if seq in unique_sequences:
			unique_sequences[seq] += 1
		else:
			for longer_seq in unique_sequences:
				##check if it is substring of sequence already in unique_sequences
				if seq in longer_seq:
					unique_sequences[longer_seq] += 1
					is_substring = True
					"substring found!"
					break
			##if still unique (i.e. not a substring), add to dictionary of counts
			if not is_substring:	
				unique_sequences[seq] = 1
	
	confident_unique_seqs = list(k for k, v in unique_sequences.items() if v > 1)
	
	print "... Total unique sequence variants: %s" % len(unique_sequences)
	print "... Unique sequence variants associated with >1 template molecules: %s" % len(confident_unique_seqs) 
	
	###STEP 6: Name and align sequences, generate fasta output
	print "[INFO 6]: Creating output fasta '%s.unique_sequences.fas' of aligned unique sequence variants with headers '>Sample_Name-#AssociatedTemplates'" % prefix
	##Header is sample prefix+number of templates associated with sequence
	unique_records = []
	for s in unique_sequences:
		sequence = Seq(s)
		name = '%s-%s'  % (prefix, unique_sequences[s])
		seq_rec = SeqRecord(sequence, id=name, description='')
		unique_records.append(seq_rec)
	
	##Add reference sequence 
	reference_rec = SeqRecord(Seq(reference_seq), id='reference-sequence', description='')
	unique_records.append(reference_rec)
	
	##Align the unique sequences and order with most templates at top
	aligned_unique_records = get_alignment(unique_records)
	aligned_unique_records.sort(reverse=True)
	
	##Write to fasta file
	SeqIO.write(aligned_unique_records, '%s.unique_sequences.fas' % prefix, 'fasta')
	
	##Wrap up
	#print "Done."	