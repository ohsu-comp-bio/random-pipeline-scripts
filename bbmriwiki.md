# BBMRI Pipeline

http://www.bbmriwiki.nl/wiki/PipelineCommands

Below is the current GoNL pipeline commands along with rough execution times when available on the cluster as of January 25th 2011:

This analysis pipeline consists of two parts. The first one is based on lane analysis. Afterwards lanes are merged and analysis per lane is conducted.

## Lane analysis 
**L1. Fastq to check read quality**

```bash
/tools/FastQC/fastqc /data/filename_1.fq.gz \ -Dfastqc.output_dir=/output \ -Dfastqc.unzip=false

/tools/FastQC/fastqc /data/filename_2.fq.gz \ -Dfastqc.output_dir=/output \ -Dfastqc.unzip=false
```
**L2. Align both pairs of FASTQ files using BWA **(~6.5h)**[[BR]] **

```bash
/tools/bwa-0.5.8c_patched/bwa aln \ /resources//hg19/indices/b37_1kg.fa \ /data/filename_1.fq.gz -t 4 \ -f /output/filename.b37_1kg.1.sai

/tools/bwa-0.5.8c_patched/bwa aln \ /resources//hg19/indices/b37_1kg.fa \ /data/filename_2.fq.gz -t 4 \ -f /output/filename.b37_1kg.2.sai
```
**L3. Create SAM file from both .sai files using BWA sampe **(~3.5h)**[[BR]] **

```bash
/tools/bwa_45_patched/bwa sampe -P -p illumina -i lane -m sample_id -l library \ /resources//hg19/indices/b37_1kg.fa \ /output/filename.b37_1kg.1.sai /output/filename.b37_1kg.2.sai \ /data/filename_1.fq.gz /data/filename_2.fq.gz \ -f /output/filename.b37_1kg.sam
```
**L4. Convert SAM to BAM file using Picard **(~0.5h)**[[BR]] **

```bash
java -jar -Xmx3g /tools/picard-tools-1.32/SamFormatConverter.jar \ INPUT=/output/filename.b37_1kg.sam \ OUTPUT=/output/filename.b37_1kg.bam \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local
```
**L5. Sort BAM file and create corresponding index **(~5h)**[[BR]] **

```bash
java -jar -Xmx3g /tools/picard-tools-1.32/SortSam.jar \ INPUT=/output/filename.b37_1kg.bam \ OUTPUT=/output/filename.b37_1kg.sorted.bam \ SORT_ORDER=coordinate \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar -Xmx3g /tools/picard-tools-1.32/BuildBamIndex.jar \ INPUT=/output/filename.b37_1kg.sorted.bam \ OUTPUT=/output/filename.b37_1kg.sorted.bam.bai \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local
```
**L6. Calculate QC metrics on alignment using Picard** (~3.5h, excluding coverage - currently not working)

This Step is currently updated as:

 * Picard !CalculateHsMetrics should be replaced by GATK !DepthOfCoverage
 * !MeanQualityByCycle and !QualityScoreDistribution should be moved after the recalibration step

```bash
java -jar -Xmx4g /tools/picard-tools-1.32/CollectAlignmentSummaryMetrics.jar \ I=/output/filename.b37_1kg.sorted.bam \ O=/output/filename.b37_1kg.AlignmentSummaryMetrics \ R=/resources//hg19/indices/b37_1kg.fa \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar /tools/picard-tools-1.32/CollectGcBiasMetrics.jar \ R=/resources//hg19/indices/b37_1kg.fa \ I=/output/filename.b37_1kg.sorted.bam \ O=/output/filename.b37_1kg.GcBiasMetrics \ CHART=/output/filename.b37_1kg.GcBiasMetrics.pdf \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar /tools/picard-tools-1.32/CollectInsertSizeMetrics.jar \ I=/output/filename.b37_1kg.sorted.bam \ O=/output/filename.b37_1kg.CollectInsertSizeMetrics \ H=/output/filename.b37_1kg.CollectInsertSizeMetrics.pdf \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar /tools/picard-tools-1.32/MeanQualityByCycle.jar \ I=/output/filename.b37_1kg.sorted.bam \ O=/output/filename.b37_1kg.MeanQualityByCycle \ CHART=/output/filename.b37_1kg.MeanQualityByCycle.pdf \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local java -jar /tools/picard-tools-1.32/QualityScoreDistribution.jar \ I=/output/filename.b37_1kg.sorted.bam \ O=/output/filename.b37_1kg.[wiki:QualityScoreDistribution] \ CHART=/output/filename.b37_1kg.[wiki:QualityScoreDistribution].pdf \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar /tools/picard-tools-1.32/BamIndexStats.jar \ INPUT=/output/filename.b37_1kg.sorted.bam \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar -Xmx3g /tools/picard-tools-1.32/CalculateHsMetricsWholeGenome.jar \ INPUT=/output/filename.b37_1kg.sorted.bam \ OUTPUT=/output/filename.b37_1kg.HsMetrics \ BAIT_INTERVALS=/resources//hg19/intervals/GoNL.interval_list \ TARGET_INTERVALS=/resources//hg19/intervals/GoNL.interval_list \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local
```
**L7. Mark duplicates after alignment** (~2.5h)

```bash
java -Xmx4g -jar /tools/picard-tools-1.32/MarkDuplicates.jar \ INPUT=/output/filename.b37_1kg.sorted.bam \ OUTPUT=/output/filename.b37_1kg.dedup.bam \ METRICS_FILE=/output/filename.b37_1kg.dedup.metrics \ REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local
```
**L8. Realignment using knowns only** (~2.5h)

```bash
java -Djava.io.tmpdir=/local -Xmx8g -jar \ /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO -T IndelRealigner \ -U ALLOW_UNINDEXED_BAM -I /output/filename.b37_1kg.dedup.bam \ -targetIntervals /resources//hg19/intervals/realign_intervals_hg19_b37_1kg.intervals \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -[B:indels,VCF B:indels,VCF] /resources//hg19/indels/1kg.pilot_release.merged.indels.sites./hg19.b37_1kg.vcf \ -o /output/filename.b37_1kg.realigned.bam -knownsOnly -LOD 0.4 -compress 0
```
**L9. Fix mates after realignment** (~2h)

```bash
java -jar -Xmx6g /tools/picard-tools-1.32/FixMateInformation.jar \ INPUT=/output/filename.b37_1kg.realigned.bam \ OUTPUT=/output/filename.b37_1kg.matefixed.bam \ SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/local
```
**L10. Calculate covariates before realignment ** (~12h - Note that this should be improved using more cores. In test at the moment)

```bash
java -jar -Xmx2g /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO \ -T CountCovariates -U ALLOW_UNINDEXED_BAM \ -R /resources//hg19/indices/b37_1kg.fa \ --DBSNP /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -I /output/filename.b37_1kg.matefixed.bam \ -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \ -recalFile /output/filename.b37_1kg.matefixed.covariate_table.csv
```
**L11. Recalibrate mate fixed and realigned alignment** (~5h)

```bash
java -jar -Xmx4g /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO \ -T TableRecalibration -U ALLOW_UNINDEXED_BAM \ -R /resources//hg19/indices/b37_1kg.fa -I /output/filename.b37_1kg.matefixed.bam \ --recal_file /output/filename.b37_1kg.matefixed.covariate_table.csv \ --out /output/filename.b37_1kg.recal.bam
```
**L12. Sort and index recalibrated alignment** (~5h)

```bash
java -jar -Xmx3g /tools/picard-tools-1.32/SortSam.jar \ INPUT=/output/filename.b37_1kg.recal.bam \ OUTPUT=/output/filename.b37_1kg.recal.sorted.bam \ SORT_ORDER=coordinate \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local

java -jar -Xmx3g /tools/picard-tools-1.32/BuildBamIndex.jar \ INPUT=/output/filename.b37_1kg.recal.sorted.bam \ OUTPUT=/output/filename.b37_1kg.recal.sorted.bam.bai \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local
```
**L13. Calculate covariates after realignment and recalibration** (~12h Note that this should be improved using more cores. In test at the moment)

```bash
java -jar -Xmx2g /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO \ -T CountCovariates -U ALLOW_UNINDEXED_BAM \ -R /resources//hg19/indices/b37_1kg.fa \ --DBSNP /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -I /output/filename.b37_1kg.recal.sorted.bam \ -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \ -recalFile /output/filename.b37_1kg.recal.covariate_table.csv
```
**L14. Analyze covariates before and after** (Currently not working)

```bash
java -jar -Xmx4g /tools/Sting/dist/AnalyzeCovariates.jar -l INFO \ -resources /resources//hg19/indices/b37_1kg.fa \ --recal_file /output/filename.b37_1kg.matefixed.covariate_table.csv \ -outputDir /output/filename.b37_1kg.recal.stats_before/ \ -Rscript ${rscript} -ignoreQ 5

java -jar -Xmx4g /tools/Sting/dist/AnalyzeCovariates.jar -l INFO \ -resources /resources//hg19/indices/b37_1kg.fa \ --recal_file /output/filename.b37_1kg.recal.covariate_table.csv \ -outputDir /output/filename.b37_1kg.recal.stats_after/ \ -Rscript ${rscript} -ignoreQ 5
```
= Sample analysis =
**S1. Merging three lanes to one** (~7h)

```bash
java -jar -Xmx4g /tools/picard-tools-1.32/MergeSamFiles.jar \ INPUT/output/filename.b37_1kg.recal.sorted.ba \ INPUT/output/filename2.b37_1kg.recal.sorted.ba \ INPUT/output/filename3.b37_1kg.recal.sorted.ba \ ASSUME_SORTED=true USE_THREADING=true \ TMP_DIR=/local MAX_RECORDS_IN_RAM=30000000 \ OUTPUT=/output/sample_id.b37_1kg.bam SORT_ORDER=coordinate \ VALIDATION_STRINGENCY=SILENT
```
== To perform QC on initial data SNP calling is performed. This procedure consists of three steps. ==
**S2. SNP calling using Unified Genotyper**

```bash
java -Xmx4g -jar /tools/Sting/dist/GenomeAnalysisTK.jar \ -l INFO -T UnifiedGenotyper -I /output/sample_id.b37_1kg.bam \ --out /output/sample_id.b37_1kg.qc_check_snps.vcf \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -stand_call_conf 30.0 -stand_emit_conf 10.0
```
**S3. Filter SNPs using variant filtration**

```bash
java -jar -Xmx4g /tools/Sting/dist/GenomeAnalysisTK.jar \ -l INFO -T VariantFiltration \ -[B:variant,VCF B:variant,VCF] /output/sample_id.b37_1kg.qc_check_snps.vcf \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ --out /output/sample_id.b37_1kg.qc_check_snps.filtered.vcf \ --maskName InDel --clusterWindowSize 10 \ --filterName GATK_standard \ --filterExpression "AB > 0.75 && DP > 40 || DP > 100 || MQ0 > 40 || SB > -0.10"||
```
**S4. Generate statistics on called SNPs using variant eval **

```bash
java -Xmx4g -jar /tools/Sting/dist/GenomeAnalysisTK.jar \ -T VariantEval -R /resources//hg19/indices/b37_1kg.fa \ -l INFO \ -[B:eval,VCF B:eval,VCF] /output/sample_id.b37_1kg.qc_check_snps.filtered.vcf \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -o /output/sample_id.b37_1kg.qc_check_snps.filtered.eval
```
== After this QC check the indel and SNP calling takes place using the next commands. ==
**S5. Create targets for realignment**

```bash
java -Xmx4g -jar -Djava.io.tmpdir=/local /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO \ -T RealignerTargetCreator \ -I /output/sample_id.b37_1kg.bam \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -[B:indels,VCF B:indels,VCF] /resources//hg19/indels/1kg.pilot_release.merged.indels.sites./hg19.b37_1kg.vcf \ -o /output/sample_id.b37_1kg.realign.intervals
```
**S6. Realignment using created targets, dbSNP and indels from 1kg project**

```bash
java -Djava.io.tmpdir=/local –Xmx6g -jar \ /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO -T IndelRealigner \ -I /output/sample_id.b37_1kg.bam \ -targetIntervals /output/sample_id.b37_1kg.realign.intervals \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -[B:indels,VCF B:indels,VCF] /resources//hg19/indels/1kg.pilot_release.merged.indels.sites./hg19.b37_1kg.vcf \ --out /output/sample_id.b37_1kg.realigned.bam \ -maxReads 500000
```
**S7.  Fix mates after realignment and create corresponding BAM index**

```bash
java -jar -Xmx6g /tools/picard-tools-1.32/FixMateInformation.jar \ INPUT=/output/sample_id.b37_1kg.realigned.bam \ OUTPUT=/output/sample_id.b37_1kg.matesfixed.bam \ SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/local

java -jar -Xmx3g /tools/picard-tools-1.32/BuildBamIndex.jar \ INPUT=/output/sample_id.b37_1kg.matesfixed.bam \ OUTPUT=/output/sample_id.b37_1kg.matesfixed.bam.bai \ VALIDATION_STRINGENCY=LENIENT \ TMP_DIR=/local
```
**S8. Indel calling on realigned BAM file**

```bash
java -Xmx4g -jar /tools/Sting/dist/GenomeAnalysisTK.jar -l INFO \ -T IndelGenotyperV2 -I /output/sample_id.b37_1kg.matesfixed.bam \ --out /output/sample_id.b37_1kg.indels.vcf \ --bedOutput /output/sample_id.b37_1kg.indels.bed \ -R /resources//hg19/indices/b37_1kg.fa \ -verbose /output/sample_id.b37_1kg.indels.verboseoutput.txt
```
**S9. Filter indels **

```bash
perl /tools/filterSingleSampleCalls.pl \ --calls /output/sample_id.b37_1kg.indels.bed > /output/sample_id.b37_1kg.indels.filtered.bed \ --max_cons_av_mm 3.0 --max_cons_nqs_av_mm 0.5 --mode ANNOTATE
```
**S10. SNP calling on realigned BAM file**

```bash
java -Xmx4g -jar /tools/Sting/dist/GenomeAnalysisTK.jar \ -l INFO -T UnifiedGenotyper -I /output/sample_id.b37_1kg.matesfixed.bam \ --out /output/sample_id.b37_1kg.snps.vcf \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -stand_call_conf 30.0 -stand_emit_conf 10.0
```
**S11. Create indel mask used in SNP calling python **

```bash
/tools/makeIndelMask.py \ /output/sample_id.b37_1kg.indels.filtered.bed 10 \ /output/sample_id.b37_1kg.indels.mask.bed
```
**S12. Filter SNP calls using the indel mask, dbsnp and indels from 1kg project**

```bash
java -jar -Xmx4g /tools/Sting/dist/GenomeAnalysisTK.jar \ -l INFO -T VariantFiltration \ -[B:variant,VCF B:variant,VCF] /output/sample_id.b37_1kg.snps.vcf \ -[B:mask,Bed B:mask,Bed] /output/sample_id.b37_1kg.indels.mask.bed \ -R /resources//hg19/indices/b37_1kg.fa \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ --out /output/sample_id.b37_1kg.snps.filtered.vcf \ --maskName InDel --clusterWindowSize 10 \ --filterName GATK_standard \ --filterExpression "AB > 0.75 && DP > 40 || DP > 100 || MQ0 > 40 || SB > -0.10"||
```
**S13. Variant eval on detected SNPs**

```bash
java -Xmx4g -jar /tools/Sting/dist/GenomeAnalysisTK.jar \ -T VariantEval -R /resources//hg19/indices/b37_1kg.fa \ -l INFO \ -[B:eval,VCF B:eval,VCF] /output/sample_id.b37_1kg.snps.filtered.vcf \ -D /resources//hg19/dbsnp/dbsnp_129_b37_b37_1kg.rod \ -o /output/sample_id.b37_1kg.snps.filtered.eval
```
