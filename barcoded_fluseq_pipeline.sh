#! /bin/bash
# ~/scripts/barcoded-flu-seq/barcoded_fluseq_pipeline.sh -i $sample_prefix

##Paths to software dependencies
#PEAR: http://sco.h-its.org/exelixis/web/software/pear/
pear_path='/Users/pclangat/Software/pear/bin/pear'
quasr_path='/Users/pclangat/scripts/scripts_from_others/QUASR/QUASR_v7.02/readsetProcessor.jar'
r_path='/usr/local/bin/R'
barcode_filter_path='/Users/pclangat/scripts/barcoded-flu-seq/primerID_pipeline.py'

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

##Create directory with prefix name if it doesn't already exist, move files into it, move into it
echo "[INFO 1]: Creating directory for sample '$prefix'."
mkdir -p $sample_path

echo "[INFO 2]: Copying sample files into sample directory."
for i in `ls ${PWD}/$prefix*`;
do
	test -f $i && echo ...$i && cp $i $sample_path/
done

echo "[INFO 3]: Checking sample directory for read files."
cd $sample_path

##Get fastq read files
read1="_R1_001.fastq.gz"
read2="_R2_001.fastq.gz"

if [ -f $prefix$read1 ]; then
	echo "... Found read1 '$prefix$read1'"
	if [ -f $prefix$read2 ]; then
		echo "... Found read2 '$prefix$read2'"
	else
		echo "[ERROR]: Read2 '$prefix$read2' not found."
	fi
else
	echo "[ERROR]: Read1 '$prefix$read1' not found."
	exit 1
fi

##Pair fastq files
echo '[INFO 4]: Pairing read files using PEAR.'
$pear_path -f $prefix$read1 -r $prefix$read2 -o $prefix

##Check if pairing worked
if [ -s $prefix.assembled.fastq ]; then
	echo '... Paired file found OK.'
else
	echo "[ERROR]: Paired output file is empty. Problem with pairing occurred."
	exit 1
fi 

##Tidy up/zipup fastq files
echo '[INFO 5]: Cleaning up and compressing fastq files'.
for fastq in ./$prefix.*.fastq
do
	if [ -s $fastq ]; then	
			gzip -f $fastq
	else
		rm $fastq
	fi
done


##Run QC on individual fastq
length=400
quality=20
echo "[INFO 6]: Performing quality control using QUASR with -l $length -m $quality options:"
java -jar $quasr_path -i $prefix.assembled.fastq.gz -o $prefix -q -l $length -m $quality -g -R $r_path 

##Check if QC worked
if [ -s $prefix.qc.fq ]; then
	echo '... QC reads found.'
else
	echo "[ERROR]: QC output file is empty/not found. Problem with QC occurred."
	exit 1
fi

##Barcode filtering
echo '[INFO 7]: Filtering reads with sequencing barcodes.'
python $barcode_filter_path $prefix

##Finish
echo '[INFO 8]: Cleaning up and compressing remaining fastq files'.
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
echo "Done."
mv log.$prefix.txt $sample_path