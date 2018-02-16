#! /bin/bash
# ~/scripts/barcoded-flu-seq/wrapper_script.sh -i <prefix_list_file>

## Takes file of prefixes as input, runs barcoded_fluseq_pipeline.sh for each prefix and outputs log file

##Help message
#if [ "$#" = "0" ]; then
#	echo "[USAGE]: barcoded_fluseq_pipeline.sh <options>
#-h   [None]		Print this help message
#-i * [String]		Input file of samples"
#	exit 1
#fi

##Accepting and checking user input
while getopts ":i:h" flag; do
	case $flag in
		h) echo "[USAGE]: assemblyPipeline_map.sh <options>
	-h   [None]		Print this help message
	-i * [String]		Input file of samples"
	    exit 1
		;;
		i) infile=$OPTARG
		;;
		:) echo "[ERROR]: Flag -$OPTARG requires an argument. See -h for details."
		   exit 1
		;;
		?) echo "[INFO]: Ignoring invalid option: -$OPTARG. See -h for details."
		   exit 1
	esac
done

if [ "$infile" = '' ]; then
	echo '[ERROR]: Incorrect usage. See -h for details.'
	exit 1
fi

## Read listfile

read1="_R1_001.fastq.gz"
read2="_R2_001.fastq.gz"
paired=".pair.fastq.gz"

while read -r line
do
	prefix="$line"
	echo "Merging paired reads for $prefix..."
	/Users/pclangat/Software/bbmap/fuse.sh in1=$prefix$read1 in2=$prefix$read2 out=$prefix$paired fusepairs pad=2 overwrite=t
	#/Users/pclangat/scripts/barcoded-flu-seq/barcoded_fluseq_pipeline.sh -i $prefix > "log.$prefix.txt"
done < "$infile"
	