# barcoded-flu-seq
Processing barcoded (Primer ID) influenza HA sequencing data

General Dependencies: Need to download the following and set your paths in the '.sh' script

##Paths to software dependencies
#PEAR: http://sco.h-its.org/exelixis/web/software/pear/
pear_path='/Users/pclangat/Software/pear/bin/pear'
quasr_path='/Users/pclangat/scripts/scripts_from_others/QUASR/QUASR_v7.02/readsetProcessor.jar'
r_path='/usr/local/bin/R'
barcode_filter_path='/Users/pclangat/scripts/barcoded-flu-seq/primerID_pipeline.py'

Python dependencies: need to install the following for the '.py' script
Requirements: 
requires pip: brew install pip
requires biopython: pip install biopython
requires regex: pip install regex

To run: 
~/scripts/barcoded-flu-seq/barcoded_fluseq_pipeline.sh -i $sample_prefix
