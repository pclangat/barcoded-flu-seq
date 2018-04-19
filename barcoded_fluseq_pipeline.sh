#! /bin/bash
# ~/scripts/barcoded-flu-seq/barcoded_fluseq_pipeline.sh -i $sample_prefix

##Paths to software dependencies
#PEAR: http://sco.h-its.org/exelixis/web/software/pear/
pear_path='/software/pathogen/external/apps/usr/bin/pear'
quasr_path='/software/pathogen/external/apps/usr/local/QUASR/readsetProcessor.jar'
r_path='/software/R-3.4.2/bin/R'
barcode_filter_path='/nfs/users/nfs_p/pl6/scripts/barcoded-flu-seq/primerID_pipeline.py'

##Help message
if [ "$#" = "0" ]; then
	echo "[USAGE]: barcoded_fluseq_pipeline.sh <options>
-h   [None]		Print this help message
-i * [String]		Input file prefix"
	exit 1
fi

##Accepting and checking user input
while getopts ":i:h" flag; do
	case $flag in
		h) echo "[USAGE]: assemblyPipeline_map.sh <options>
	-h   [None]		Print this help message
	-i * [String]		Input file prefix"
	    exit 1
		;;
		i) prefix=$OPTARG
		;;
		:) echo "[ERROR]: Flag -$OPTARG requires an argument. See -h for details."
		   exit 1
		;;
		?) echo "[INFO]: Ignoring invalid option: -$OPTARG. See -h for details."
		   exit 1
	esac
done

if [ "$prefix" = '' ]; then
	echo '[ERROR]: Incorrect usage. See -h for details.'
	exit 1
fi

##Sample path
sample_path="${PWD}/$prefix/"

echo "[INFO]: Checking sample directory for read files."
cd $sample_path

##Get fastq read files
reads=".fastq"

if [ -f *$reads ]; then
	echo -e "\t... Found basecalled reads."
	
else
	echo "[ERROR]: Reads not found."
	exit 1
fi

##Initial quality check
echo "[INFO 4]: Initial quality check of reads using QUASR."
#java -jar $quasr_path -i $prefix$read1 -o $prefix$r1 -g -w 300 -R $r_path
#java -jar $quasr_path -i $prefix$read2 -o $prefix$r2 -g -w 300 -R $r_path 
echo "...currently skipping this step."

##Unzip fastq files
echo "[INFO 5a]: Gunzipping fastq read files."
gunzip $prefix$read1
gunzip $prefix$read2

##Pair fastq files
echo "[INFO 5b]: Pairing read files using PEAR."
$pear_path -f $prefix$read1unzip -r $prefix$read2unzip -o $prefix

##Check if pairing worked
if [ -s $prefix.assembled.fastq ]; then
#if [ -s $prefix.pair.fastq.gz ]; then
	echo -e "\t... Paired file found OK."
else
	echo "[ERROR]: Paired output file is empty. Problem with pairing occurred."
	exit 1
fi 

##Run QC on paired fastq
##Set min length, phred quality, and expected paired reads length
length=300
quality=20
pairedlength=600
echo "[INFO 6]: Performing quality control using QUASR with -l $length -m $quality options:"
#java -jar $quasr_path -i $prefix.assembled.fastq -o $prefix -q -l $length -m $quality -z -g -w $length -R $r_path
#java -jar $quasr_path -i $prefix.pair.fastq.gz -o $prefix -q -l $length -m $quality -g -w $pairedlength -R $r_path  
#java -jar $quasr_path -i $prefix.assembled.fastq -o $prefix -q -l $length -m $quality -g -w $pairedlength -R $r_path
java -jar $quasr_path -i $prefix.assembled.fastq -o $prefix -q -l $length -m $quality
 
##Check if QC worked
#if [ -s $prefix.qc.fq.gz ]; then
if [ -s $prefix.qc.fq ]; then
	echo -e '\t... QC reads found.'
else
	echo "[ERROR]: QC output file is empty/not found. Problem with QC occurred."
	exit 1
fi

##Tidy up/zipup fastq files
echo '[INFO 7]: Cleaning up and compressing fastq files'.
for fastq in ./*.fastq
do
	if [ -s $fastq ]; then	
			gzip -f $fastq
	else
		rm $fastq
	fi
done

##Barcode filtering
echo '[INFO 8]: Running PrimerID error correction/template counting on reads with sequencing barcodes.'
python $barcode_filter_path $prefix

##Finish
echo '[INFO 9]: Cleaning up and compressing remaining fastq files'.
for fastq in ./*.fq
do
	if [ -s $fastq ]; then	
			gzip -f $fastq
	else
		rm $fastq
	fi
done

mkdir ./intermediate_fastqs
for gz in ./*.gz
do
	mv $gz ./intermediate_fastqs
done

cd ..
mv log.$prefix.o $sample_path
mv log.$prefix.e $sample_path
echo "Done."