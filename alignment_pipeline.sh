#!/bin/bash

### Specify number of worker threads.
DTHREADS=20
CTHREADS=4
MAX_READS=14000000

### Program directories
BIN_DIR="/opt/installed"
PICARD="$BIN_DIR/picard/picard-tools-1.110"
#GATK="$BIN_DIR/GATK/GenomeAnalysisTK-3.1.jar"
GATK="/mnt/lustre1/CompBio/bin/GATK-3.5.jar"

### Project directories
REF_GENOME_DIR="/mnt/lustre1/CompBio/genomic_resources/genomes/rat/"
PROJ="/mnt/lustre1/users/balter/horton_chondrocyte_chipseq_balter/alignment"
if [ ! -e $PROJ ]; then
    mkdir $PROJ
fi
RAW="$PROJ/raw_bam"
LOGS="$PROJ/logs"

SAMPLE=$1

if [ $SAMPLE == "99999" ]; then
    rm -rf $PROJ/{bwa, sortsam, markdups, rtc, ir, br, pr, hc, temp}
#    rm -rf $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/ir $PROJ/br $PROJ/pr $PROJ/hc $PROJ/temp
    rm $LOGS/*
    exit
fi

### Create directory structure.
mkdir -p $PROJ/{bwa, sortsam, markdups, rtc, ir, br, pr, hc, temp}
#mkdir -p $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/ir $PROJ/br $PROJ/pr $PROJ/hc $PROJ/temp

### Set read groups
READGROUPS="@RG\tID:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE"

### Copy resource files to temporary sample directories.  This helps prevent multiple processes from accessing the same files and causing problems.
### Probably just need to provide these files on each node, especially if pipelines will be used by general public.

### Can probably use condor transfer files to automatically move them and then delete them.
#######################
mkdir -p $PROJ/temp/$SAMPLE
cp $MILLS $PROJ/temp/$SAMPLE
cp $MILLS.idx $PROJ/temp/$SAMPLE
cp $DBSNP $PROJ/temp/$SAMPLE
cp $DBSNP.idx $PROJ/temp/$SAMPLE
cp $TG $PROJ/temp/$SAMPLE
cp $TG.idx $PROJ/temp/$SAMPLE
#######################

### Align with BWA
if [ ! -e $PROJ/sortsam/$SAMPLE\_sortsam.bam.md5 ]; then
    bwa mem 
        -t $DTHREADS 
        -k 19 
        -w 100 
        -d 100 
        -r 1.5
        -c 10000 
        -A 1 
        -B 4 
        -O 6 
        -E 1 
        -L 5 
        -U 17 
        -T 30 
        -M 
        -R $READGROUPS 
            $HG $RAW/$SAMPLE\_R1.fastq $RAW/$SAMPLE\_R2.fastq > $PROJ/bwa/$SAMPLE.bam

### Sort alignments
    java -jar -Xms36g -Xmx54g $PICARD/SortSam.jar
        SO=coordinate
        I=$PROJ/bwa/$SAMPLE.bam
        O=$PROJ/sortsam/$SAMPLE\_sortsam.bam
        MAX_RECORDS_IN_RAM=$MAX_READS 
        CREATE_MD5_FILE=true
fi

### Mark Duplicates
if [ ! -e $PROJ/markdups/$SAMPLE\_markdups.bam.md5 ]; then
    java -jar -Xms36g -Xmx54g $PICARD/MarkDuplicates.jar 
        MAX_FILE_HANDLES=200 
        MAX_RECORDS_IN_RAM=$MAX_READS 
        I=$PROJ/sortsam/$SAMPLE\_sortsam.bam
        O=$PROJ/markdups/$SAMPLE\_markdups.bam
        M=$PROJ/markdups/$SAMPLE.metrics
        CREATE_INDEX=true
        CREATE_MD5_FILE=true
fi

if [ ! -e $PROJ/ir/$SAMPLE\_ir.bam.md5 ]; then

### Realign Target Creator
    java -jar $GATK
        -R $HG
        -T RealignerTargetCreator
        -I $PROJ/markdups/$SAMPLE\_markdups.bam 
        -nt $DTHREADS
        -o $PROJ/rtc/$SAMPLE.interval_list 
        -known:indels,vcf $PROJ/temp/$SAMPLE/mills.vcf 
        --disable_auto_index_creation_and_locking_when_reading_rods

### Indel Realigner
    java -jar $GATK
        -R $HG 
        -T IndelRealigner 
        -I $PROJ/markdups/$SAMPLE\_markdups.bam 
        -targetIntervals $PROJ/rtc/$SAMPLE.interval_list 
        -maxInMemory $MAX_READS 
        -known:indels,vcf $PROJ/temp/$SAMPLE/mills.vcf 
        -o $PROJ/ir/$SAMPLE\_ir.bam 
        --generate_md5 
        --disable_auto_index_creation_and_locking_when_reading_rods
fi
