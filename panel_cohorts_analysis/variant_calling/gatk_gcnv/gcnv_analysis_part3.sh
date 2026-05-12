#!/bin/bash

##Job settings
#PBS -V
#PBS -N p3_HVOL

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=8:00:00
#PBS -J 1-1805

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=HVOL
#########################

WORKDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/${cohort}_temp
CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
REFDICT=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.dict
TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/genelists/${cohort}_genes.bed
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/${cohort}

SAMPLE_NO=$((${PBS_ARRAY_INDEX} - 1))

singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" PostprocessGermlineCNVCalls \
--model-shard-path $WORKDIR/variant-calling-model \
--calls-shard-path $WORKDIR/variant-calling-calls \
--allosomal-contig chrX \
--allosomal-contig chrY \
--contig-ploidy-calls $WORKDIR/ploidy-calls \
--sample-index ${SAMPLE_NO} \
--output-genotyped-intervals $WORKDIR/sample${SAMPLE_NO}.intervals.vcf.gz \
--output-genotyped-segments $WORKDIR/sample${SAMPLE_NO}.segments.vcf.gz \
--output-denoised-copy-ratios $WORKDIR/sample${SAMPLE_NO}.copyratios.vcf.gz \
--sequence-dictionary $REFDICT

#keep only variant sites in the genes of interest
#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/bcftools

bcftools view \
--regions-file $TARGET \
--include 'GT!="0/0" & GT!="0" & GT!="./."' \
$WORKDIR/sample${SAMPLE_NO}.segments.vcf.gz \
--output-type z \
--output-file $OUTDIR/sample${SAMPLE_NO}.FINAL.vcf.gz

#remove output files if there are no variants of interest
VARIANT_COUNT=$(zcat $OUTDIR/sample${SAMPLE_NO}.FINAL.vcf.gz | grep -v '^#' | wc -l)

if [ $VARIANT_COUNT -eq 0 ]
then
    rm -f $OUTDIR/sample${SAMPLE_NO}*
fi