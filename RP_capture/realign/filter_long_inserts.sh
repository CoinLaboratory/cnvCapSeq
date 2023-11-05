#!/bin/sh

#This step requires SAMtools and Bowtie2

chrom=$1
SAMPLE=$2
INPUT=$3
OUTPUT=$4
JARDIR=$5
start=$6
stop=$7
insertThresh=$8
AUX=$9

file=$SAMPLE

mkdir -p $OUTPUT/split_samples
echo "Extracting abnormal read pairs..."
java -Xmx3800m -jar $JARDIR/SplitByInsertSize.jar $SAMPLE $chrom $INPUT $OUTPUT $start $stop $insertThresh

mkdir -p $OUTPUT/reads
echo "Creating fastq files for realignment..."
java -Xmx3800m -jar $JARDIR/SamToFastq.jar I=$OUTPUT/split_samples/$file.abnormal.bam F=$OUTPUT/reads/$file.abnormal._1.fq F2=$OUTPUT/reads/$file.abnormal._2.fq

mkdir -p $OUTPUT/realignments
echo "Realigning using Bowtie2..."
bowtie2 -t -x $AUX/chr$chrom.fa -1 $OUTPUT/reads/$file.abnormal._1.fq -2 $OUTPUT/reads/$file.abnormal._2.fq --rg-id "$SAMPLE\_sorted.bam" --rg "LB:LB_$SAMPLE\_sorted.bam" --rg "SM:SM_$SAMPLE" --rg "PL:illumina" --rg "PU:PU_$SAMPLE" --rg "DS:DS_$SAMPLE" --rg "CN:CN_$SAMPLE" -a --very-sensitive --maxins 300000 -p 4 | samtools view -uhS - | samtools sort - $OUTPUT/realignments/$file.abnormal.realign.bowtie

samtools index $OUTPUT/realignments/$file.abnormal.realign.bowtie.bam

mkdir -p $OUTPUT/unique_realignments
echo "Selecting optimal realignment..."
java -Xmx3800m -jar $JARDIR/BestAltMapping.jar $SAMPLE $chrom $OUTPUT 150 300 $insertThresh
samtools index $OUTPUT/unique_realignments/$file.abnormal.realign.bowtie.unique.bam

mkdir -p $OUTPUT/final_filtered
java -Xmx3800m -jar $JARDIR/FilterAbnormalInsertSizes.jar $SAMPLE $chrom $OUTPUT

mkdir -p $OUTPUT/processed
echo "Merging and sorting realignments..."
samtools merge -u - $OUTPUT/split_samples/$SAMPLE.concordant.bam $OUTPUT/final_filtered/$SAMPLE.abnormal.final.filtered.bam | samtools sort - $OUTPUT/processed/$SAMPLE.filtered

samtools index $OUTPUT/processed/$SAMPLE.filtered.bam
