#!/bin/bash

##Job settings
#PBS -V
#PBS -N gcnv_p1_HVOL

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=2:00:00
#PBS -J 1-1805

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=HVOL
#########################

#DATADIR=$EPHEMERAL/TEST/$cohort
DATADIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/$cohort
CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
REFGENOME=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
BEDFILE=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.bed

#preprocess intervals from bed if not already done

if [ ! -e /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.interval_list ]
then
    singularity exec $CONTAINERDIR/gatk_latest.sif \
    gatk PreprocessIntervals \
    -R $REFGENOME \
    -L $BEDFILE \
    --padding 250 \
    --bin-length 0 \
    -imr OVERLAPPING_ONLY \
    -O /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.interval_list
fi

#create filelists for analysis and define PBS array index
if [ ! -e /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt ]
then
    ls $DATADIR/*.bam > /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt
fi
#ls $DATADIR/*.bam > /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/TESTfiles$cohort.txt

BAMFILE=$(sed -n ${PBS_ARRAY_INDEX}p /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt)
#BAMFILE=$(sed -n ${PBS_ARRAY_INDEX}p /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/TESTfiles$cohort.txt)

BAMFILENAME=$(basename $BAMFILE)
SAMPLEID="${BAMFILENAME%.bam}" 

#collect read count info
WORKDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/${cohort}_temp
mkdir -p $WORKDIR

singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk CollectReadCounts \
-L /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.interval_list \
-R $REFGENOME \
-imr OVERLAPPING_ONLY \
-I $BAMFILE \
--format TSV \
-O $WORKDIR/$SAMPLEID.tsv

#annotate intervals
if [ ! -e /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.annotated.tsv ]
then
    singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
    gatk AnnotateIntervals \
    -L /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.interval_list \
    -R $REFGENOME \
    -imr OVERLAPPING_ONLY \
    -O /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.annotated.tsv
fi