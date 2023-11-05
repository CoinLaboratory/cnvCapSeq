#!/bin/sh

#This step requires SAMtools

export JARDIR=$1
export OUTPUTROOT=$2
export INPUT=$3
export OUTPUT=$OUTPUTROOT/Processed
export TMPDIR=$OUTPUTROOT/Temp

sample=$4
maxInsertSize=$5
start=$6
stop=$7

mkdir -p $OUTPUT
mkdir -p $TMPDIR

#First we sort and index the input BAM file using SAMtools

samtools sort $INPUT/$sample.filtered.bam $TMPDIR/$sample.sorted
samtools index $TMPDIR/$sample.sorted.bam

#Then we split the BAM file into its individual libraries and normalize the libraries
#The "ProcessInsertSize" class takes 4 arguments:
#-1 the path of the input BAM file
#-2 the path of the temporary directory where the intermediate files will be stored
#-3 the output path
#-4 the maximum allowed RP distance (defaults to 300000bp)
#-5 the start coordinate
#-6 the end coordinate

java -Xms3800m -Xmx3800m -cp $JARDIR/ProcessInsertSize.jar paired_end.ProcessInsertSize $TMPDIR/$sample.sorted.bam $TMPDIR $OUTPUT $maxInsertSize $start $stop

#Delete the intermediate files
rm -rf $TMPDIR/full_split
rm -rf $TMPDIR/sorted_ranked_signed_libraries
rm -rf $TMPDIR/normalized_to_each_other
rm -rf $TMPDIR/final_libraries
