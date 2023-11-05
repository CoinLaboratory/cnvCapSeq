#!/bin/sh

export JARDIR=$1
export OUTPUTROOT=$2
export INPUT=$OUTPUTROOT/Consolidated/final_adjusted_inserts/
export OUTPUT=$OUTPUTROOT/Final_Merged/
export chr=$3
export start=$4
export stop=$5
export window=$6

mkdir -p $OUTPUT

#This class merges and samples the individual processed
#read pair files into the format recognized by cnvCapSeq.
#The larger the sample set the more memory it requires.

#This step can also make use of a toExclude.txt file in the ./RP_capture/RP_pipeline
#directory, in order to exclude regions that are known to have failed the capture
#step. The toExclude.txt file is not mandatory and can be empty.

#The "CombineAndSampleInsertCap" class takes 9 arguments:
#-1 the path of the folder that contains the individually processed files
#-2 the path of the output folder
#-3 the size of the sampling window (defaults to 20, increase if you have memory problems)
#-4 the chromosome that we are processing
#-5 the coordinate from which to start (defaults to 1)
#-6 the maximum coordinate that will be processed.
#-7 the path of the file that contains possible regions to be excluded (absent here)
#-8 specifies the tolerance around the excluded region breakpoints (defaults to 50)
#-9 the start coordinate of the region of interest if different from the previously defined
#-10 the minimum percentage of read pairs (relative to the region average) that need to span each position

java -Xmx3800m -cp $JARDIR/CombineAndSampleInsertCap.jar cap.CombineAndSampleInsertCap $INPUT $OUTPUT $window $chr $start $stop $JARDIR/toExclude.txt 50 $start 0.3
