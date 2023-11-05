#!/bin/sh

#RD_pipeline_capture script.

#First we need to define the working directory and the input directory.
#The input directory should only contain the BAM (and possibly the BAI index) files.
#The working directory is were all the intermediate and final results will be saved
#We assume each BAM file contains reads from one sample, for one chromosome.

chr=$1 #chromosome to perform analysis on
start=${2:-0} #start coordinate
stop=$3 #end coordinate
INPUT=$4 #input folder (containing the BAM files)

export OUTPUTDIR=./pre_processing_results/RD_results
export JARDIR=./RD_capture/
export INPUT=${INPUT:="./bam/"}
export FINAL=./CNV_detection_capture/data/RD

mkdir -p $OUTPUTDIR
mkdir -p $FINAL

#First, we process the BAM files by index, sorting and piling them up.

for f in $INPUT/*.bam
do
	echo "Processing file: $(basename $f) ..."	
	sh $JARDIR/process_one_capture.sh $OUTPUTDIR $INPUT $f
done


#The next step is to sample (according to the windowSize parameter)
#and smoothe the read depth using the singular value decomposition.

windowSize=50

echo "Smoothing Read Depth..."
sh $JARDIR/capture_smoothe_all.sh $OUTPUTDIR $JARDIR $chr $start $stop $windowSize


#The final step is to merge the sampled, smoothed files.

echo "Merging files..."
sh $JARDIR/merge_all_capture.sh $JARDIR $OUTPUTDIR $chr $start $stop $windowSize

cp $OUTPUTDIR/Final_Merged/chr$chr"_"$start"_"$stop/*.zip $FINAL/$chr.zip

echo "Results copied to ./CNV_detection/data/RD folder"

