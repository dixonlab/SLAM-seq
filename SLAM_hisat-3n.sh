#!/bin/sh

FASTQ_DIR=~/data1/2022_06_17 # path to fastq.gz files
DATA_DIR=~/SLAM/RAD21-dTAG/timecourse/novaseq # output directory

CPUS=8

# dir/basename used in hisat-3n-build
GEN_DIR=~/RPE1-SNPs/HISAT-3N_SLAMseq/hs38d5_RPE1_Nsub.with-decoy.TC

hisat=/pbld/raidix/home/tpopay/apps/hisat-3n # path to hisat-3n (see https://daehwankimlab.github.io/hisat2/hisat-3n/)
JAVA2="java -Xmx60g -jar"
PICARD=~/apps/picard-tools-2.1.1/picard.jar

mkdir $DATA_DIR/tmp
TMP=$DATA_DIR/tmp

for NAME in TP470_S23 TP471_S24 # list of fastq file names excluding extension
do
READ1=$NAME\_R1_001
READ2=$NAME\_R2_001

cd $hisat
OUT=$DATA_DIR/hisat-3n-FR/$NAME
mkdir -p $OUT
mkdir -p $OUT/map
mkdir -p $OUT/counts

./hisat-3n --base-change T,C -5 10 -x $GEN_DIR --rna-strandness FR --unique-only --directional-mapping --threads $CPUS -q -1 $FASTQ_DIR/$READ2\.fastq.gz -2 $FASTQ_DIR/$READ1\.fastq.gz -S $OUT/map/hisat2_aligned.sam
samtools sort $OUT/map/hisat2_aligned.sam -o $OUT/map/hisat2_aligned.sorted.bam -O bam
$JAVA2 $PICARD MarkDuplicates INPUT=$OUT/map/hisat2_aligned.sorted.bam OUTPUT=$OUT/map/hisat2_aligned.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=$OUT/map/metrics.$NAME.txt TMP_DIR=$TMP

samtools index $OUT/map/hisat2_aligned.nodup.bam
samtools sort -n -o $OUT/map/hisat2_aligned.nodup.name_sorted.bam $OUT/map/hisat2_aligned.nodup.bam

done

