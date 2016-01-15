**STILL BEING REFORMATTED FROM SOURCE WIKI MARKUP**

# Utah Genome Project
(From http://weatherby.genetics.utah.edu/UGP/wiki/index.php/UGP_Variant_Pipeline_1.2.1 Jan 14 2016)

 Jan. 2015 
 Variant Calling Pipeline Version 1.2.1

## Software Versions
*[http://srynobio.github.io/cApTUrE/ cApTUrE] is a lightweight NGS pipeline, created for the  Utah Genome Project (UGP)
*BWA: 0.7.10
*Picard: 1.127 (Broad version)
*GATK: 3.3-0
*SamTools: 0.2.0
*FastQC v0.10.1

## Data Source

Data sets used for the variant calling pipeline come from the Broad GSA (GATK) group as the 'GATK resource bundle version 2.8

```
wget -r ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/b37
```

### Reference Genome (GRCh37):
* human_g1k_v37_decoy.fasta

### Call region file generated from NCBI
* GRCh37 GFF3

### VCF files for RealignerTargetCreator knowns and dbsnp for BaseRecalibrator.
* known_indel: Mills_and_1000G_gold_standard.indels.b37.vcf
* known_indel: 1000G_phase1.indels.b37.vcf
* known_dbsnp: dbsnp_137.b37.vcf

## Background Files
* We have created 1000Genomes (BWA mem/GATK 3.0+) background files to be ran concurrently with the GenotypeGVCFs step.

Groups Currently completed:
*CEU
*GBR
*FIN
''Version 1.0.5 background files have been updated to show only the indviduals of each group, not the file names.''

''This is a complete list of the background individuals for run completed > 1.0.5 [http://weatherby.genetics.utah.edu/UGP/wiki/index.php/Background-individuals]''

 If you would like a copy of the current files, we have made a public AWS s3 bucket
 Using [http://s3tools.org/s3cmd s3cmd] execute the following command: 
 '''s3cmd get s3://ugp-1k-backgrounds --recursive'''

 Alternatively to access the files without s3cmd the following use the following URLs:
 http://s3-us-west-2.amazonaws.com/ugp-1k-backgrounds/CEU_mergeGvcf.vcf
 http://s3-us-west-2.amazonaws.com/ugp-1k-backgrounds/CEU_mergeGvcf.vcf.idx
 http://s3-us-west-2.amazonaws.com/ugp-1k-backgrounds/FIN_mergeGvcf.vcf
 http://s3-us-west-2.amazonaws.com/ugp-1k-backgrounds/FIN_mergeGvcf.vcf.idx
 http://s3-us-west-2.amazonaws.com/ugp-1k-backgrounds/GBR_mergeGvcf.vcf
 http://s3-us-west-2.amazonaws.com/ugp-1k-backgrounds/GBR_mergeGvcf.vcf.idx

Resource files for VariantRecalibrator_SNP
*hapmap_3.3.b37.vcf
*1000G_omni2.5.b37.vcf
*1000G_phase1.snps.high_confidence.b37.vcf

Resource files for VariantRecalibrator_INDEL
*Mills_and_1000G_gold_standard.indels.b37.vcf
*1000G_phase1.indels.b37.vcf

== Sequencing ==

This pipeline is designed for 100 bp (or greater) Illumina HiSeq PE exome or WGS sequence data with Sanger/Illumina 1.9 quality encoding, and uses Illumina naming convention found here [http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm]

=== Validate File Integrity with md5sum ===

An md5sum signature should be provided for each FastQ file by the sequencing center.  After the file has been downloaded, locally check the md5sum to be sure that no data corruption occurred during the file transfer.

 md5sum file.fastq > file_local.md5
 diff file_local.md5 file_provided.md5
 Now the pipeline runs md5_check to validation the results and will quit if errors are found. 

If the md5sum signature differs from that provided for the file:
* Check to be sure you have the correct file.
* Check if the md5sum was calculated on that compressed or uncompressed file by the provider and be sure to do the same with the local copy.
* Try the download again.
* Contact the sequence provider.

== FastQ File Analyses ==

 fastqc Sample1_L1_R1.txt

From the sumamry.txt report we check 
* FAIL sections

From the fastqc_data.txt file we check the following values:
* Encoding (must be Sanger / Illumina 1.9)
* Total Sequences (Currently set to 30000000)
* Filtered Sequences (Currently set to less then 5)
* Sequence length (must be >= 100 bp)
* %GC (45 < x < 55)
* Total Duplicate Percentage (Currently set to 60.0)
''The pipeline now runs fastqc_check and output these result into QC-report.txt.''

== Indexing ==

The following indexing is required using BWA, Picard and SamTools.  GATK requires all three.  However this step only needs to be done once "per-machine".

*BWA
 bwa index -a bwtsw human_g1k_v37_decoy.fasta

*Picard
 java -jar CreateSequenceDictionary.jar R=human_g1k_v37_decoy.fasta O=human_g1k_v37_decoy.dic

*SamTools
 samtools faidx human_g1k_v37_decoy.fasta

== Alignment ==

Align reads to the genome with bwa.

The 'BWA-mem' program will find the reference coordinates of the input reads (independent of their mate-pair). The following parameters are those used by the 1KG project and GATK for aligning Illumina data.

 bwa mem -M -R "read group" human_g1k_v37_decoy.fasta Sample1_L1_R1.fq Sample1_L1_R2.fq | samtools view -bSho BAM_FILE -

== BAM File Analyses and Processing ==

Alignment BAM files are improved in various ways to help increase the quality and speed of subsequent variant calling steps.

=== BAM Sorting ===

 java -jar -XX:ParallelGCThreads=# -Xmx#g -Djava.io.tmpdir=/tmp /usr/local/picard/dist/picard.jar SortSam 
    INPUT=*.bam 
    OUTPUT=*_sorted.bam 
    VALIDATION_STRINGENCY=LENIENT 
    SORT_ORDER=coordinate 
    MAX_RECORDS_IN_RAM=5000000 
    CREATE_INDEX=True 
    &> SortSam.log-#

=== Merge lane level BAMs to individual === 

* This step only needs to run if you have multiple lanes per sample.

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp /usr/local/picard/dist/picard.jar MergeSamFiles.jar
     INPUT=#*_sorted.bam
     INPUT=#*_sorted.bam
     INPUT=[ ... ]
     OUTPUT=*.bam                          
     VALIDATION_STRINGENCY: LENIENT
     MAX_RECORDS_IN_RAM: 5000000
     CREATE_INDEX: True
     SORT_ORDER: coordinate
     ASSUME_SORTED: True 
     USE_THREADING: True
     &> MergeSamFiles.log-#

=== Mark Duplicates ===

Remove PCR/Optical duplicate reads

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp /usr/local/picard/dist/picard.jar MarkDuplicates 
    INPUT=*_sorted.bam 
    OUTPUT=*_sorted_Dedup.bam 
    METRICS_FILE=*_sorted_Dedup.metrics 
    VALIDATION_STRINGENCY=LENIENT 
    MAX_RECORDS_IN_RAM=5000000 
    ASSUME_SORTED=True 
    CREATE_INDEX=True  
    &> MarkDuplicates.log-#

<span style="background:#FFFF00">Currently all duplicates are flagged to allow GATK to handle them.</span>

=== BAM Quality Control ===

''At this point the pipeline has broken the tasks into chromosomal regions.''

*CollectMultipleMetrics
 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp /usr/local/picard/dist/picard.jar CollectMultipleMetrics
    INPUT=*_sorted_Dedup_realign_chr#_recal.bam
    OUTPUT=*_sorted_Dedup_realign_chr#_recal.metrics 
    VALIDATION_STRINGENCY=LENIENT 
    PROGRAM=QualityScoreDistribution  
    REFERENCE_SEQUENCE=human_g1k_v37_decoy.fasta 
    &> CollectMultipleMetrics.log-#

*idxstats
 samtools idxstats [dedup bam files] > dedup-bamfile.stats

''Now the pipeline will take idxstats ouput and check for unmapped reads.''

=== Local Realignment of Indels ===

*RealignerTargetCreator

''This is where the tasks are broken into chromosomal regions.''

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T RealignerTargetCreator 
    -R human_g1k_v37_decoy.fasta 
    -I *_sorted_Dedup.bam 
    --num_threads #  
    --known Mills_and_1000G_gold_standard.indels.b37.vcf
    --known 1000G_phase1.indels.b37.vcf  
    -L chr#_region_file.list 
    -o *_chr#_realign.intervals 
    &> RealignerTargetCreator.log-#

*IndelRealigner

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T IndelRealigner 
    -R human_g1k_v37_decoy.fasta 
    -I *_sorted_Dedup.bam 
    -L chr#_region_file.list
    -targetIntervals *_chr#_realign.intervals 
    -known Mills_and_1000G_gold_standard.indels.b37.vcf 
    -known 1000G_phase1.indels.b37.vcf 
    -o *_sorted_Dedup_realign_chr#.bam 
    &> IndelRealigner.log-1

=== BaseRecalibration & PrintReads ===

*BaseRecalibrator

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T BaseRecalibrator 
    -R /human_g1k_v37_decoy.fasta 
    -I *_sorted_Dedup_realign_chr#.bam 
    --num_cpu_threads_per_data_thread #  
    --knownSites dbsnp_137.b37.vcf
    --knownSites Mills_and_1000G_gold_standard.indels.b37.vcf 
    --knownSites 1000G_phase1.indels.b37.vcf  
    -o *_sorted_Dedup_realign_chr#_recal_data.table 
    &> BaseRecalibrator.log-#

*PrintReads

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T PrintReads 
    -R human_g1k_v37_decoy.fasta 
    -I *_sorted_Dedup_realign_chr#.bam
    --num_cpu_threads_per_data_thread #  
    -BQSR *_sorted_Dedup_realign_chr#_recal_data.table 
    -o *_sorted_Dedup_realign_chr#_recal.bam 
    &> PrintReads.log-#
 
<span style="background:#FFFF00">

== Variant Calling ==

=== HaplotypeCaller ===

*Now HaplotypeCaller handels SNP and INDEL calls.

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T HaplotypeCaller 
    -R human_g1k_v37_decoy.fasta 
    --min_base_quality_score 20 
    --variant_index_parameter 128000 
    --emitRefConfidence GVCF 
    --standard_min_confidence_threshold_for_calling 30.0 
    --num_cpu_threads_per_data_thread # 
    --variant_index_type LINEAR 
    --standard_min_confidence_threshold_for_emitting 30.0  
    -I *_sorted_Dedup_realign_chr*_recal.bam 
    -L chr#_region_file.list
    -o chr#_region_*.raw.snps.indels.gvcf 
    &> HaplotypeCaller.log-#

=== CatVariants ===

''Collect all individual gvcf files.''

 java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants 
    -R human_g1k_v37_decoy.fasta 
    --assumeSorted  
    -V chr#_region_*.raw.snps.indels.vcf 
    -V [ ... ]
    -out *.raw.snps.indels.gCat.vcf 
    &> CatVariants.log-*

=== CombineGVCFs ===

''Collect all the chromosomal cat files.''
''This is also the step where mergeGvcf are collected for future runs.''

 java -jar -Xmx#g -XX:ParallelGCThreads=# GenomeAnalysisTK.jar  
    -T CombineGVCFs 
    -R human_g1k_v37_decoy.fasta 
    --variant *.raw.snps.indels.gCat.vcf
    --variant [ ... ]
    -o cApTUrE_*_final_mergeGvcf.vcf 
    &> CombineGVCF.log-*

=== GenotypeGVCFs ===

 java -jar -Xmx#g -XX:ParallelGCThreads=# GenomeAnalysisTK.jar 
    -T GenotypeGVCFs 
    -R human_g1k_v37_decoy.fasta 
    --num_threads #  
    --variant cApTUrE_*_final_mergeGvcf.vcf 
    --variant CEU_mergeGvcf.vcf 
    --variant GBR_mergeGvcf.vcf 
    --variant FIN_mergeGvcf.vcf 
    -L chr#_region_file.list 
    -o cApTUrE_*_#_genotyped.vcf 
    &> GenotypeGVCF.log-#

=== Combine_Genotyped ===

''Collect all the individual genotype steps.''

 java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants 
    -R human_g1k_v37_decoy.fasta 
    --assumeSorted  
    -V cApTUrE_*_#_genotyped.vcf 
    -V [ ... ] 
    -out cApTUrE_*_genotyped.vcf 
    &> Combine_Genotyped.log-#

=== VariantRecalibrator ===

*SNP Recalibration

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=tmp GenomeAnalysisTK.jar  
    -T VariantRecalibrator 
    -R human_g1k_v37_decoy.fasta 
    --minNumBadVariants 5000 
    --num_threads #  
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf 
    -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf 
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf 
    -an QD 
    -an MQRankSum 
    -an ReadPosRankSum 
    -an FS 
    -an InbreedingCoeff 
    -input cApTUrE_*_genotyped.vcf 
    -recalFile cApTUrE_*_snp_recal 
    -tranchesFile cApTUrE_*_snp_tranches 
    -rscriptFile cApTUrE_*_snp_plots.R 
    -mode SNP 
    &> VariantRecalibrator_SNP.log-#

*INDEL Recalibration

 java -jar -Xmx#g -XX:ParallelGCThreads=# -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T VariantRecalibrator 
    -R human_g1k_v37_decoy.fasta 
    --minNumBadVariants 5000 
    --num_threads #  
    -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf 
    -resource:1000G,known=false,training=true,truth=true,prior=10.0 1000G_phase1.indels.b37.vcf 
    -an MQRankSum 
    -an ReadPosRankSum 
    -an FS 
    -input cApTUrE_*_genotyped.vcf 
    -recalFile cApTUrE_*_indel_recal 
    -tranchesFile cApTUrE_*_indel_tranches 
    -rscriptFile cApTUrE_*_indel_plots.R 
    -mode INDEL 
    &> VariantRecalibrator_INDEL.log-#

=== ApplyRecalibration ===

*SNP Apply
 
 java -jar -Xmx#g -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T ApplyRecalibration 
    -R human_g1k_v37_decoy.fasta 
    --ts_filter_level 99.5 
    --excludeFiltered 
    --num_threads #  
    -input cApTUrE_*_genotyped.vcf 
    -recalFile cApTUrE_*_snp_recal 
    -tranchesFile cApTUrE_*_snp_tranches 
    -mode SNP 
    -o cApTUrE_*_recal_SNP.vcf 
    &> ApplyRecalibration_SNP.log-*

*INDEL Apply

 java -jar -Xmx#g -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar
    -T ApplyRecalibration 
    -R human_g1k_v37_decoy.fasta 
    --ts_filter_level 99.0 
    --excludeFiltered 
    --num_threads #  
    -input cApTUrE_*_genotyped.vcf 
    -recalFile cApTUrE_*_indel_recal 
    -tranchesFile cApTUrE_*_indel_tranches 
    -mode INDEL 
    -o cApTUrE_*_recal_INDEL.vcf 
    &> ApplyRecalibration_INDEL.log-*

== CombineVarients ==

These command will combine INDEL and SNP files into a single VCF file.
*CombineVarients

 java -jar -Xmx#g -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T CombineVariants 
    -R human_g1k_v37_decoy.fasta 
    --num_threads # 
    --genotypemergeoption UNSORTED 
    --variant cApTUrE_*_recal_SNP.vcf  
    --variant cApTUrE_*_recal_INDEL.vcf  
    -o cApTUrE_*_Combined.vcf 
    &> CombineVariants.log-*

*SelectVariants

 java -jar -Xmx#g -Djava.io.tmpdir=/tmp GenomeAnalysisTK.jar 
    -T SelectVariants 
    -R human_g1k_v37_decoy.fasta 
    --variant cApTUrE_*_Combined.vcf  
    -select "DP > 100" 
    -o cApTUrE_*_Final+Backgrounds.vcf 
    &> SelectVariants.log-*

== Variant File QC ==

<span style="background:#FFFF00">
Quality Metrics on variants
*Ti/Tv Ratio (2.1 for WGS ~2.8 for exome)
*HapMap concordance
*SNV/Indel Counts
*Rare variant enrichment
*DP
*Q
*GQ
