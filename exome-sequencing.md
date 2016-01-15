# Exome Sequencing
From https://www.biostars.org/p/1268/ Jan 14 2016

Here is our basic pipeline. We are novices, but this seems to work for us. Please investigate all of the tools carefully, as (I'll repeat) I am not an expert in this analysis. Note that I am using mostly default values for the analysis; advice I gotten from experienced hands is it's routine to use the intersection of multiple SNP calling methods, and that indel calling is still an art. Note that this is set up for single reads, not paired reads, but the same basic pipeline should apply.

## PREPROCESS:
1. Index human genome (Picard), we used HG18 from UCSC.
1. Convert Illumina reads to Fastq format
1. Convert Illumina 1.6 read quality scores to standard Sanger scores

## FOR EACH SAMPLE:
1. Align samples to genome (BWA), generates SAI files.
1. Convert SAI to SAM (BWA)
1. Convert SAM to BAM binary format (SAM Tools)
1. Sort BAM (SAM Tools)
1. Index BAM (SAM Tools)
1. Identify target regions for realignment (Genome Analysis Toolkit)
1. Realign BAM to get better Indel calling (Genome Analysis Toolkit)
1. Reindex the realigned BAM (SAM Tools)
1. Call Indels (Genome Analysis Toolkit)
1. Call SNPs (Genome Analysis Toolkit)
1. View aligned reads in BAM/BAI (Integrated Genome Viewer)

## Sample Script:

```bash
/bin/BWA/bwa aln -f /output/FOO.sai -t 3 /seq/REFERENCE/human_18.fasta /seq/FQ/FOO.sanger.fq

/bin/BWA/bwa samse -f /output/FOO.sam /seq/REFERENCE/human_18.fasta /output/FOO.sai /seq/FQ/FOO.sanger.fq

/bin/SAMTOOLS/samtools import /seq/REFERENCE/human_18.fasta /output/FOO.sam /output/FOO.bam

/bin/SAMTOOLS/samtools sort /output/FOO.bam /output/FOO.sorted

/bin/SAMTOOLS/samtools index /output/FOO.sorted.bam /output/FOO.sorted.bam.bai

java -jar /bin/GTK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /seq/REFERENCE/human_18.fasta -I /output/FOO.sorted.bam  -o /output/FOO.intervals

java -jar /bin/GTK/GenomeAnalysisTK.jar -T IndelRealigner -R /seq/REFERENCE/human_18.fasta -I /output/FOO.sorted.bam -targetIntervals /output/FOO.intervals --output /output/FOO.sorted.realigned.bam

/bin/SAMTOOLS/samtools index /output/FOO.sorted.realigned.bam /output/FOO.sorted.realigned.bam.bai

java -jar /bin/GTK/GenomeAnalysisTK.jar -T IndelGenotyperV2 -R /seq/REFERENCE/human_18.fasta -I /output/FOO.sorted.realigned.bam -O /output/FOO_indel.txt --verbose -o /output/FOO_indel_statistics.txt

java -jar /bin/GTK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /seq/REFERENCE/human_18.fasta -I /output/FOO.sorted.realigned.bam -varout /output/FOO.geli.calls -vf GELI -stand_call_conf 30.0 -stand_emit_conf 10.0 -pl SOLEXA
```
