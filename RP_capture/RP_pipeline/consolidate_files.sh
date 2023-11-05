#!/bin/sh

export JARDIR=$1
export OUTPUTROOT=$2
export OUTPUT=$OUTPUTROOT/Consolidated
export start=$3
export stop=$4

mkdir -p $OUTPUT

#The next step is to put the libraries of the BAM file back together.
#The "ConsolidateFiles" class takes 4 arguments:
#-1 the path of the directory that contains the processed library files
#-2 the output path
#-3 the start coordinate
#-4 the end coordinate

java -Xmx3800m -server -cp $JARDIR/ConsolidateFiles.jar paired_end.ConsolidateFiles $OUTPUTROOT/Processed $OUTPUT $start $stop
