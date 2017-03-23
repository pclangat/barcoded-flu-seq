#!/usr/bin/python
## ~/scripts/primerID_pipeline.py <sample_name>

## Takes as input the sample name --> finds fwd and reverse reads in fastq.gz file
## STEP ONE: Pairs forward and reverse. Requires PEAR.
## STEP TWO: Does QC via QUASR. Requires QUASR.

##Import modules
import sys, os, regex
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

###METHODS

##Ensure correct directory created and all files moved
#def ensure_correct_path(p):
#	current_path = os.getcwd()
#	sample_path = ("%s/%s" % (current_path, p))
#	
#	##Create sample directory
#	if not os.path.exists(sample_dir):
#		os.makedirs(sample_dir)
#		print("Creating directory for sample %s" % p)
#	
#	##Move files to directory

def get_sequence_file(prfx, sufx, type):
	if 'fastq' in type:
		ext = '.fastq.gz'
	elif 'fasta' in type:
		ext = '.fas'
	else:
		print("Error finding sequence file type: %s" % type)
		sys.exit(1)
	file_name = prfx + sufx + ext
	return file_name

def run_PEAR_pairing(prfx, f1, f2):
	print("Attempting to pair %s and %s" % (f1, f2))
	pair_command = pear_path + ' -f ' + f1 + ' -r ' + f2 + ' -o ' + prfx
	os.system(pair_command)
	return

def run_QUASR_qc(fastq, prfx, length, phred):
	print("Attempting to run QC on " + fastq)
	qc_command = 'java -jar ' + quasr_path+'/QUASR_v7.02/readsetProcessor.jar -i ' + fastq + ' -o ' + prfx + ' -q -l ' + str(length) + ' -m ' + str(phred) + ' -g -R ' + r_path
	os.system(qc_command)
	return

def check_rev_primer_match(r_primer, seq):
	match = 'none'
	pos = seq.find(r_primer)
	
	## if match found
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
		#else:
			##if still no match, try searching for partial match
				
	return pos, seq, match
		
def fuzzy_substring(needle, haystack):
	"""Calculates the fuzzy match of needle in haystack,
	using a modified version of the Levenshtein distance
	algorithm.
	The function is modified from the levenshtein function
	in the bktree module by Adam Hupp"""
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
		print row1
	print min(row1)
	print (0.3*m)
	if min(row1) < (0.3*m):
		primer_pos = row1.index(min(row1)) - m
	else:
		primer_pos = -1
	return primer_pos
				
###INITIALISATIONS
##Paths to software
quasr_path = "/Users/pl6/scripts/scripts_from_others/QUASR/"
r_path = "/usr/bin/R"
pear_path = "/Users/pl6/Software_downloaded/pear/bin/pear"

f_sufx = "_R1_001"
r_sufx = "_R2_001"

##These will become user-accepted later
min_length = 125
min_phred = 20

fwd_primer = "CGGGGAAAATATGCAACAATCCT"
rev_primer = "GAGGGTTTCACTTGGACTGGG"
full_rev_primer = "GAGGGTTTCACTTGGACTGGGNNNNANNNNANNNNAAAGCAGTGGTATCAACGCA"
pattern = rev_primer+'[A-Z]{4}A[A-Z]{4}A[A-Z]{4}A'
pattern = regex.compile(pattern)

###MAIN
if __name__ == '__main__':
	##Parse user input
	try:
		prfx = sys.argv[1]
	except:
		print('[USAGE]: %s sample_prefix' % sys.argv[0] )
		sys.exit(1)

	##Get fastq files
	#fastq1 = get_sequence_file(prfx, f_sufx, 'fastq')
	#fastq2 = get_sequence_file(prfx, r_sufx, 'fastq')
	
	##Pair the fastq files
	#run_PEAR_pairing(prfx, fastq1, fastq2)
	#paired_fastq = prfx + '.assembled.fastq'
	
	##Run QC on paired fastq
	#run_QUASR_qc(paired_fastq, prfx, min_length, min_phred)
	qc_paired_fastq = prfx + '.qc.fq'
	qc_paired_fasta = prfx + '.fas'
	
	##Convert qc fastq files to fasta files
	#SeqIO.convert(qc_paired_fastq, "fastq", qc_paired_fasta, "fasta")

	##Load qc FASTAs as dicts. SeqIO only reads up to the space, so now the names should match in 1 and 2
	read_records = SeqIO.index(qc_paired_fasta, 'fasta')
	reads_dict = {}
	
	##Only keeps what is needed from SeqIO record in dict
	for read in read_records:
		reads_dict[read] = read_records[read].seq
	read_records.close()
	
	##Find reverse primer (gene portion) in read sequence, count whether full or partial or no match
	matches = []
	barcodes=0
	barcoded_seqs = defaultdict(dict)
	for read in reads_dict:
		sequence = reads_dict[read]
		primer_pos, oriented_sequence, match_depth = check_rev_primer_match(rev_primer, sequence)
		matches.append(match_depth)
		reads_dict[read] = oriented_sequence
		if primer_pos > 0:
			id = oriented_sequence[primer_pos:]
			id = regex.findall(pattern, str(id))
			if id:	
				barcodes+=1
				barcode = id[0].split(rev_primer)[1]
				if oriented_sequence in barcoded_seqs[barcode]:
					barcoded_seqs[barcode][oriented_sequence] += 1 
				else:
					barcoded_seqs[barcode][oriented_sequence] = 1 
	print barcoded_seqs
	print "Full matches to reverse primer: %s" % matches.count('full')
	print "Partial matches to reverse primer: %s" % matches.count('partial')
	print "Reverse primer not found: %s" % matches.count('none')
	print "Intact barcodes found: %s" % barcodes
	
	##Barcode filtering
	
	##Wrap up
	print "Done."	
		
		
		
