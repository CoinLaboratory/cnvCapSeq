#!/bin/sh

#RP_pipeline script for use with a PBS queuing system.

#First, we need to define the working directory and the input directory. The 
#input directory should only contain the BAM files. BAI index files will created as
#necessary. BAM files are assumed to contain reads from one sample, for one chromosome.
#The working directory is were all the intermediate and final results will be saved

chr=$1 #chromosome to perform analysis on
start=${2:-1} #start coordinate
stop=$3 #end coordinate
INPUT=$4 #input folder (containing the BAM files)

export OUTPUTDIR=./pre_processing_results/RP_results
export JARDIR=./RP_capture/
export INPUT=${INPUT:="./bam/"}
export FINAL=./CNV_detection_capture/data/RP
export AUX=./aux

mkdir -p $OUTPUTDIR

#This step extracts abnormal read pairs from the bam file
#and performs sensitive local realignment in an effort to
#salvage them. The results are re-incorporated into the original
#bam file that will is used in downstream processing.

sh $JARDIR/realign/realign.sh $chr $start $stop $INPUT $OUTPUTDIR $JARDIR $AUX

export OUTPUT=$OUTPUTDIR
export JARDIR=./RP_capture/RP_pipeline
export INPUT=$OUTPUTDIR/realign/processed


#Then, we split the BAM files into their original sequencing libraries
#which are then quantile normalised to a Gaussian N(200,15).
#This step requires SAMtools (http://samtools.sourceforge.net/)

maxInsertSize=300000
for f in $INPUT/*.bam
do
	echo "Normalizing $(basename $f .filtered.bam)"	
	sh $JARDIR/process_one.sh $JARDIR $OUTPUT $INPUT $(basename $f .filtered.bam) $maxInsertSize $start $stop
done

#The next step is to consolidate the normalized library files.

echo "Consolidating files..."
sh $JARDIR/consolidate_files.sh $JARDIR $OUTPUT $start $stop 


#Finally, we merge and sample the normalised files. This final command should
#be executed when the previous job has finished. The user can specify which
#genomic region will be included in the merged file, and how dense the sampling
#will be (in bp), using the window variable. The default is 20bp.

windowSize=100
echo "Merging files..."
sh $JARDIR/merge_all.sh $JARDIR $OUTPUT $chr $start $stop $windowSize

cp $OUTPUTDIR/Final_Merged/chr$chr"_"$start"_"$stop/*.zip $FINAL/$chr.zip

echo "Results copied to ./CNV_detection/data/RP folder"
