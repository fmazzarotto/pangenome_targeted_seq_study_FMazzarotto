#!/bin/bash

##Job settings
#PBS -V
#PBS -N p2_HG002

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=2:00:00
#PBS -J 1-169

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=SARCWES_forHG002
#########################

echo "Part 2 of GATK gCNV analysis of HG002 and 167 control WES (+HG005) samples has started..."
echo $(date)

#set max n of cpus for JAVA/GATK to use, to avoid PBS killing of the job
export GATK_JAVA_OPTS="-XX:ActiveProcessorCount=16" 

#set paths
CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
TEMPFILESDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKgCNV/tmp
REFGENOME=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta

#set array-based logic
BAMFILE=$(sed -n ${PBS_ARRAY_INDEX}p /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt)
BAMFILENAME=$(basename $BAMFILE)
SAMPLEID="${BAMFILENAME%.bam}" 

#collect read counts (array-based on all 169 samples)
singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CollectReadCounts \
-L $TEMPFILESDIR/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid_padded_GATKgCNVanalysis.interval_list \
-R $REFGENOME \
-imr OVERLAPPING_ONLY \
-I $BAMFILE \
--format TSV \
-O $TEMPFILESDIR/$SAMPLEID.tsv