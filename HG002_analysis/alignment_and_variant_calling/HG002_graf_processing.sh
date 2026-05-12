#!/bin/bash

##Job settings
#PBS -V
#PBS -N graf_HG002

##Job configuration
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -lwalltime=6:00:00

echo "GRAF analysis of HG002 has started..."
echo $(date)

#set paths
GRAFCONTAINER=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images/gwgs-u_1.2.sif
HG002_FOLDER=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/GIAB_SRA_WES_samples
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GRAF
REFGEN=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/pangenomes/B38.fa
DCM_VARS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/pangenomes/DCM.SV.B38.vcf.gz
TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/exome-intersect_v5covered_highconfHG002.bed

mkdir -p $OUTDIR

#realign sample with the gral aligner
singularity exec --bind /rds:/rds $GRAFCONTAINER aligner \
--vcf $DCM_VARS \
--reference $REFGEN \
--markdup \
--tmp $EPHEMERAL/graf_sort.XXXXXX \
-q $HG002_FOLDER/HG002_SRR2962669_1.fastq.gz \
-Q $HG002_FOLDER/HG002_SRR2962669_2.fastq.gz \
--read_group_platform ILLUMINA \
--read_group_sample HG002 \
--sort \
--nthreads 64 \
--index \
-o $HG002_FOLDER/HG002_GRAL.bam

#call variants with graf
singularity exec --bind /rds:/rds $GRAFCONTAINER rasm \
--graph-vcf $DCM_VARS \
--fasta $REFGEN \
--bam $HG002_FOLDER/HG002_GRAL.bam \
--vcf $OUTDIR/HG002_graf_raw.vcf \
--interval-file $TARGET \
--annotation-list all \
--total-read-limit 1000000

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

# manually apply GRAF's pipeline standard filter
conda activate bcftools
bcftools filter \
--threads 16 \
--exclude "AD_Ratio[1] < 0.20 || MBQ[1] < 15 || QD < 1 || MQRankSum < -12.5 || FS > 60" \
--mode x $OUTDIR/HG002_graf_raw.vcf \
| bcftools norm -m -both -O z -o $OUTDIR/HG002_graf_FINAL.vcf.gz

bcftools view \
--threads 8 \
-f PASS \
--include 'strlen(REF) != strlen(ALT) & (strlen(REF) >= 20 | strlen(ALT) >= 20)' $OUTDIR/HG002_graf_FINAL.vcf.gz \
| bcftools sort \
-O z \
-o $OUTDIR/HG002_graf_min20bp_FINAL.vcf.gz

bcftools view \
--threads 8 \
-f PASS \
--include '(strlen(REF) < 20 & strlen(ALT) < 20)' $OUTDIR/HG002_graf_FINAL.vcf.gz \
| bcftools sort \
-O z \
-o $OUTDIR/HG002_graf_max19bp_FINAL.vcf.gz

tabix -p vcf $OUTDIR/HG002_graf_min20bp_FINAL.vcf.gz
tabix -p vcf $OUTDIR/HG002_graf_max19bp_FINAL.vcf.gz
tabix -p vcf $OUTDIR/HG002_graf_FINAL.vcf.gz

echo "Analysis of HG002 finished!!!"
echo $(date)