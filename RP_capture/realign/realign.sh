#!/bin/sh

#The realignment step requires Bowtie2, along with a reference fasta
#file for the chromosome of interest to be placed in the ./aux directory
#The fasta file also needs to be indexed with Bowtie2.

chr=$1 #chromosome to perform analysis on
start=${2:-1} #start coordinate
stop=$3 #end coordinate
INPUT=$4 #input folder (containing the BAM files)

export OUTPUT=$5/realign
export JARDIR=$6/realign
export AUX=$7

mkdir -p $OUTPUT

insertThreshold=1000

#The insertThreshold defines the minimum insert size to be considered abnormal. 
#Read pairs with a longer insert size will be locally realigned with Bowtie2
#and the best resulting alignment will be re-incorporated into the original
#bam file, substituting the "abnormal" alignment.

for f in $INPUT/*.bam
do
	echo "Processing file: $(basename $f .bam) ..."
	sh $JARDIR/filter_long_inserts.sh $chr $(basename $f .bam) $INPUT $OUTPUT $JARDIR $start $stop $insertThreshold $AUX
done
