#!/bin/bash

##Job settings
#PBS -V
#PBS -N p4_HG002

##Job configuration
#PBS -l select=1:ncpus=32:mem=32gb
#PBS -lwalltime=8:00:00
#PBS -J 1-21

#SET COHORT TO ANALYSE
#########################
cohort=HG002    #NOTE: cohort is called SARCWES_forHG002 in part 1 script - OK to have them different
#########################

echo "Part 4 of GATK gCNV analysis of HG002 and 167 control WES (+HG005) samples has started..."
echo $(date)

#set max n of cpus for JAVA/GATK to use, to avoid PBS killing of the job
export GATK_JAVA_OPTS="-XX:ActiveProcessorCount=32" 

CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
TEMPOUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKgCNV/tmp
PLOIDYPRIORS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gatkgcnv/contig_ploidy_prior_hg38.tsv
BEDFILELIST=$TEMPOUTDIR/HG002_GATKgCNV_bedfiles_scatter.txt
ANNOT_INTERVALS=$TEMPOUTDIR/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid_padded_GATKgCNVanalysis.annotated_intervals.txt
REFGENOME=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
REFDICT=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.dict

#set array-based logic
BEDFILE=$(sed -n ${PBS_ARRAY_INDEX}p $BEDFILELIST)

#create appropriate -I flag content for the following commands
INPUTSTRING=$(find "$TEMPOUTDIR" -maxdepth 1 -name "*.tsv" -type f -print0 | xargs -0 printf '-I %s ')

#call variants
singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GermlineCNVCaller \
--run-mode COHORT \
-L $BEDFILE \
$INPUTSTRING \
--contig-ploidy-calls $TEMPOUTDIR/ploidy-calls \
--annotated-intervals $ANNOT_INTERVALS \
--interval-merging-rule OVERLAPPING_ONLY \
--output $TEMPOUTDIR \
--output-prefix variant-calling-shard${PBS_ARRAY_INDEX} \
--verbosity DEBUG