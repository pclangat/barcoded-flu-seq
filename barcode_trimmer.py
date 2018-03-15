#!/usr/bin/python
## ~/scripts/barcode_trimmer.py <sample_name>

## accepts fastq and trims barcode region, outputs barcode-trimmed 

from Bio import SeqIO
import sys, re

def get_bcprimer(rev_primer):
	## Make reverse primer shorter to account for trimming
	rev_primer = rev_primer[-10:]
	bc_pattern = 'T[A-Z]{4}T[A-Z]{4}T[A-Z]{4}'
	
	pattern = rev_primer+bc_pattern
	print("Searching for barcoded primer sequences: %s" % pattern)
	return pattern

def trim_bcprimers(records, adaptor):
	trim_count = 0
	for record in records:
		seq = str(record.seq)
		match_result = re.search(adaptor, seq)
		
		if not match_result:
			#adaptor not found, try reverse complement
			rc_record = record.reverse_complement()
			rc_record.id = record.id
			rc_record.description = record.description
			seq = str(rc_record.seq)
			match_result = re.search(adaptor, seq)
			if not match_result:
				print "..."
				yield record
			else:
				#trim off adaptor from rc
				adaptor_end = match_result.end()
				og_len = len(rc_record.seq)
				trim_len = len(rc_record[adaptor_end:].seq)
				trim_count+=1
				print("%s->%s" % (og_len, trim_len))
				yield rc_record[adaptor_end:]
		else:
			#trim off adaptor
			adaptor_end = match_result.end()
			og_len = len(record.seq)
			trim_len = len(record[adaptor_end:].seq)
			trim_count+=1
			print("%s->%s" % (og_len, trim_len))
			yield record[adaptor_end:]
	print("Trimmed %i reads." % trim_count)		

trim_count = 0	
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
	
	##STEP 2: Load BC primer sequences
	uni_x = 'TGCGTTGATACCACTGCTT'
	uni_y = 'TGTCCAGCACGCTTCAGGC'
	primer_pattern = get_bcprimer(uni_x)
	primer_pattern2 = get_bcprimer(uni_y)
	
	##STEP 3: Trim the adaptors off the reads
	trimmed_reads = trim_bcprimers(original_reads, primer_pattern)
	double_trimmed_reads = trim_bcprimers(trimmed_reads, primer_pattern2)
	
	##STEP 4: Write trimmed reads to fastq file
	count = SeqIO.write(double_trimmed_reads, '%s.barcode_trimmed.qc.fq' % prefix, 'fastq')
	print("Saved %i reads." % count)
