#!/usr/bin/python
## ~/scripts/barcode_trimmer.py <sample_name>

## accepts fastq and trims barcode region, outputs barcode-trimmed 

from Bio import SeqIO
import sys, re

def get_primer():
	#uniX: rev_primer = 'TGCGTTGATACCACTGCTT'
	#uniY: rev_primer = 'TGTCCAGCACGCTTCAGGC'
	rev_primer = 'TGTCCAGCACGCTTCAGGC'
	## Make reverse primer shorter to account for trimming
	rev_primer = rev_primer[-10:]
	bc_pattern = 'T[A-Z]{4}T[A-Z]{4}T[A-Z]{4}'
	
	pattern = rev_primer+bc_pattern
	print("Searching for barcoded primer sequences: %s" % pattern)
	return pattern

def trim_bcprimers(records, adaptor):
	len_adaptor = len(adaptor)
	for record in records:
		seq = str(record.seq)
		match_result = re.search(adaptor, seq)
		
		if not match_result:
			#adaptor not found, try reverse complement
			rc_record = record.reverse_complement()
			seq = str(rc_record.seq)
			match_result = re.search(adaptor, seq)
			if not match_result:
				yield record
			else:
				#trim off adaptor from rc
				adaptor_end = match_result.end()
				yield rc_record[adaptor_end:]
		else:
			#trim off adaptor
			adaptor_end = match_result.end()
			yield record[adaptor_end:]		
	
###MAIN
if __name__ == '__main__':
	##STEP 1: Parse user input
	try:
		prefix = sys.argv[1]
	except:
		print('[USAGE]: %s sample_prefix' % sys.argv[0] )
		sys.exit(1)

	print "[INFO]: Loading sample reads..."
	original_reads = SeqIO.parse('%s.qc.fq' % prefix, 'fastq')
	
	##STEP 2: Load BC primer sequence
	primer_pattern = get_primer()
	
	trimmed_reads = trim_bcprimers(original_reads, primer_pattern)
	count = SeqIO.write(trimmed_reads, '%s.barcode_trimmed.qc.fq', 'fastq')
	print("Saved %i reads." % count)
