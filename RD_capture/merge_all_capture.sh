#!/bin/sh

export JARDIR=$1
export OUTPUTROOT=$2
export SMOOTHED=$OUTPUTROOT/Smoothed/
export OUTPUT=$OUTPUTROOT/Final_Merged/

export chr=$3
export start=$4
export stop=$5
export window=$6

tolerance=50

mkdir -p $OUTPUT

#This class merges the smoothed, sampled files and converts
#them into the required zip format. This step can also make use of
#a toExclude.txt file in the ./RD_capture directory, in order to
#exclude regions that are known to have failed the capture
#step. The toExclude.txt file is not mandatory and can be empty.

#The "CombineCapture" class takes 9 arguments:
#-1 the path of the folder that contains the smoothed, sampled files
#-2 the path of the output folder
#-3 the size of the sampling window (defaults to 100bp)
#-4 the chromosome that we are processing
#-5 the coordinate from which to start
#-6 the maximum coordinate that will be processed
#-7 the path of the file that contains regions to be excluded (can be empty)
#-8 the tolerance around the excluded region breakpoints (in bp, defaults to 50)
#-9 the start coordinate of the region of interest if different from the previously defined

java -Xmx3800m -jar $JARDIR/CombineCapture.jar $SMOOTHED $OUTPUT $window $chr $start $stop $JARDIR/toExclude.txt $tolerance $start

echo "Done pre-processing Read Depth!"
