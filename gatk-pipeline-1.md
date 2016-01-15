# GATK Pipeline

From https://www.biostars.org/p/8237/ Jan 14 2016

## Overview
GATK is certainly the best variant caller and it's the difference is important because with indels in particular, they can be hugely different. 

I've spoken to a few people who have had issues with GATK so I've actually written a blog post with a workflow that works for GATK. The major issues seemed to be using GATK's best practices, particularly their add or replace read groups module which sometimes causes it to fall over later. By adding the read group (you can make it whatever you like) when you run bwa, then this workflow actually works. 

This was written for microbes so ploidy is set to 1 for HaplotypeCaller, but this can be set to whatever you like (do check the documentation). Originally posted here(https://approachedinthelimit.wordpress.com/2015/10/09/variant-calling-with-gatk/)


## Dependencies:

You’ll need to install picardtools, GATK, bwa and optionally, vcffilter for this workflow. Picardtools and GATK are simply .jar files so that’s no problem while you probably already have bwa installed, otherwise installation is well documented!

## The workflow

This workflow begins with short read (fastq) files and a fasta reference. First a sequence dictionary is created, short reads are aligned to the reference and read group information provided, resulting sequence alignment map (sam) file sorted and converted to binary alignment map (bam) format, duplicates marked, bam file sorted, indel targets identified, indels realigned and variants called. Simple!

For simplicity an example set of commands are provided here, where the reference is reference.fasta and the short reads are R1.fastq.gz and R2.fastq.gz. You will need to enter the paths and versions of the software being used at each step and your actual file names.

```bash
#Create sequence dictionary
java -jar~/bin/picard-tools-1.8.5/CreateSequenceDictionary.jar REFERENCE=reference.fasta OUTPUT=reference.dict

#Align reads and assign read group
bwa mem -R “@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:test\tSM:PA01” reference.fasta R1.fastq.gz R2.fastq.gz > aln.sam

#Sort sam file
java -jar ~/bin/picard-tools-1.8.5/SortSam.jar I=aln.sam O=sorted.bam SORT_ORDER=coordinate

#Mark duplicates
java -jar ~/bin/picard-tools-version/MarkDuplicates.jar I=sorted.bam O=dedup.bam METRICS_FILE=metrics.txt

#Sort bam file
java -jar ~/bin/picard-tools-version/BuildBamIndex.jar INPUT=dedup.bam

#Create realignment targets
java -jar ~/bin/GATK3.3/GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I dedup.bam -o targetintervals.list

#Indel realignment
java -jar ~/bin/GATK3.3/GenomeAnalysisTK.jar -T IndelRealigner -R PA01.fasta -I dedup.bam -targetIntervals targetintervals.list -o realigned.bam

#Call variants (HaplotypeCaller)
java -jar ~/bin/GATK3.3/GenomeAnalysisTK.jar -T HaplotypeCaller -R reference.fasta -I realigned.bam -ploidy 1 -stand_call_conf 30 -stand_emit_conf 10 -o raw.vcf

The resulting vcf file will contain your variant calls!

Then you can optionally filter the variants:

#Filter variants
~/bin/vcflib/bin/vcffilter -f ‘DP > 9’ -f ‘QUAL > 10’ raw.vcf > filtered.vcf

#Or first split the raw.vcf file into SNPs and indels:

#Extract SNPs
java -jar ~/bin/GATK3.3/GenomeAnalysisTK.jar -T SelectVariants -R reference.fasta -V raw.vcf -selectType SNP -o snps.vcf

#Extract Indels
java -jar ~/bin/GATK/GenomeAnalysisTK.jar -T SelectVariants -R reference.fasta -V raw.vcf -selectType INDEL -o indels.vcf

```
