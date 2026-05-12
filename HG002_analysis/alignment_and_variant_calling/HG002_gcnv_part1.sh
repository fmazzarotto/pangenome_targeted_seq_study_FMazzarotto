#!/bin/bash

##Job settings
#PBS -V
#PBS -N p1_HG002

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=1:00:00

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=SARCWES_forHG002
#########################

echo "Part 1 of GATK gCNV analysis of HG002 and 167 control WES (+HG005) samples has started..."
echo $(date)

#set max n of cpus for JAVA/GATK to use, to avoid PBS killing of the job
export GATK_JAVA_OPTS="-XX:ActiveProcessorCount=16" 

#set paths
CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
HG002_FOLDER=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/GIAB_SRA_WES_samples
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKgCNV
TEMPFILESDIR=$OUTDIR/tmp
REFGENOME=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid.bed

mkdir -p $TEMPFILESDIR

#preprocess intervals from bed and create file with list of samples
singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" PreprocessIntervals \
-R $REFGENOME \
-L $TARGET \
--padding 250 \
--bin-length 0 \
-imr OVERLAPPING_ONLY \
-O $TEMPFILESDIR/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid_padded_GATKgCNVanalysis.interval_list

#create file with list of samples to be analyzed
ls /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/wes_ref_controls/sarcoidosis/*.bam > /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt
readlink -f $HG002_FOLDER/HG002.bam >> /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt