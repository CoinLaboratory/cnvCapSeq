#!/bin/sh

export OUTPUTROOT=$1
export PILEUPDIR=$OUTPUTROOT/Pileup
export OUTPUT=$OUTPUTROOT/Smoothed
export JARDIR=$2
export chr=$3
export start=$4
export stop=$5
export window=$6

readDepthIndex=3
coordinateIndex=1
tolerance=50

svdToRemove=1

mkdir -p $OUTPUT

#This class samples the read depth files according to the windowSize variable
#and then smoothed them using the singular value decomposition. The larger
#the sample set and the region of interest, the more memory it requires.
#This step can also make use of a toExclude.txt file in the ./RD_capture
#directory, in order to exclude regions that are known to have failed the capture
#step. The toExclude.txt file is not mandatory and can be empty.

#The "CombineCapture" class takes 11 arguments:
#-1 the path of the folder that contains the pileup files (created by samtools)
#-2 the column index of the read depth in the pileup file
#-3 the column index of the genomic coordinate in the pileup file
#-4 the path of the file that contains regions to be excluded (because they failed the capture)
#-5 the chromosome that we are processing
#-6 the coordinate from which to start
#-7 the maximum coordinate that will be processed
#-8 the tolerance around the excluded region breakpoints (in bp, defaults to 50)
#-9 the path of the output directory
#-10 the size of the sampling window (defaults to 100bp)
#-11 how many singular values to remove (the more the smoother), recommended: 1

java -Xmx3800m -jar $JARDIR/CaptureSmoothing.jar $PILEUPDIR $readDepthIndex $coordinateIndex $JARDIR/toExclude.txt $chr $start $stop $tolerance $OUTPUT $window $svdToRemove

echo "Done smoothing."
