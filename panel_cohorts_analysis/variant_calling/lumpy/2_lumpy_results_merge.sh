#!/bin/bash

##Job settings
#PBS -V
#PBS -N lpy_mrg_DCM

##Job configuration
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -lwalltime=00:30:00

##Output log configuration
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=DCM
#########################

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

#activate conda env
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/bcftools

#set paths
VCFFILES_PATH=$EPHEMERAL/$cohort
GENE_COORDS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/genelists/${cohort}_genes.bed
OUTDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/lumpy_results

mkdir -p $OUTDIR

#concatenate VCFs
bcftools merge --missing-to-ref -O z -R $GENE_COORDS --threads 4 -o $OUTDIR/${cohort}_lumpy_raw.vcf.gz $VCFFILES_PATH/*.vcf.gz