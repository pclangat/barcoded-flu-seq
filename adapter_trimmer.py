#!/usr/bin/python
## ~/scripts/barcoded-flu-seq/adapter_trimmer.py <sample_name>

## accepts fastq and trims barcode region and WGS uni13 region, outputs adapter/barcode-trimmed 

from Bio import SeqIO
import sys, re

def get_bcprimer(rev_primer):
	## Make reverse primer shorter to account for trimming
	rev_primer = rev_primer[-10:]
	bc_pattern = 'T[A-Z]{4}T[A-Z]{4}T[A-Z]{4}'
	
	pattern = rev_primer+bc_pattern
	print("Searching for barcoded primer sequences: %s" % pattern)
	return pattern

def trim_bcprimers(records, adapter):
	trim_count = 0
	for record in records:
		seq = str(record.seq)
		match_result = re.search(adapter, seq)
		
		if not match_result:
			#adapter not found, try reverse complement
			rc_record = record.reverse_complement()
			rc_record.id = record.id
			rc_record.description = record.description
			seq = str(rc_record.seq)
			match_result = re.search(adapter, seq)
			if not match_result:
				print "..."
				yield record
			else:
				#trim off adapter from rc
				adapter_end = match_result.end()
				og_len = len(rc_record.seq)
				trim_len = len(rc_record[adapter_end:].seq)
				trim_count+=1
				print("%s->%s" % (og_len, trim_len))
				yield rc_record[adapter_end:]
		else:
			#trim off adapter
			adapter_end = match_result.end()
			og_len = len(record.seq)
			trim_len = len(record[adapter_end:].seq)
			trim_count+=1
			print("%s->%s" % (og_len, trim_len))
			yield record[adapter_end:]
	print("Trimmed %i reads." % trim_count)		

def trim_fwdprimer(records, adapter):
	trim_count = 0
	for record in records:
		seq = str(record.seq)
		match_result = re.search(adapter, seq)
		
		if not match_result:
			#second adapter not found; don't try reverse complement?
			print "..."
			yield record
		else:
			#trim off adapter
			adapter_end = match_result.end()
			og_len = len(record.seq)
			trim_len = len(record[:adapter_end].seq)
			trim_count+=1
			print("%s->%s" % (og_len, trim_len))
			yield record[:adapter_end]
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
	original_reads = SeqIO.parse('%s.fastq' % prefix, 'fastq')
	
	##STEP 2: Load BC primer sequences
	uni_x = 'TGCGTTGATACCACTGCTT'
	uni_13 = 'CCTTGTTTCTACTG'
	primer_pattern = get_bcprimer(uni_x)
	primer_pattern2 = uni_13
	
	##STEP 3: Trim the adapters off the reads
	trimmed_reads = trim_bcprimers(original_reads, primer_pattern)
	double_trimmed_reads = trim_fwdprimer(trimmed_reads, primer_pattern2)
	
	##STEP 4: Write trimmed reads to fastq file
	count = SeqIO.write(double_trimmed_reads, '%s.adapter_trimmed.fastq' % prefix, 'fastq')
	print("Saved %i reads." % count)
