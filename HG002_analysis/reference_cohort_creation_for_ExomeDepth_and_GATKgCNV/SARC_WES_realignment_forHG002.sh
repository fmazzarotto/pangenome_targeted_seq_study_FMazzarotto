#!/bin/bash

##Job settings
#PBS -V
#PBS -N WES_SARC_realign

##Job configuration
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -lwalltime=16:00:00
#PBS -J 1-167

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate alignment

#SET COHORT TO ANALYSE
#########################
cohort=SARCWES_forHG002
#########################

#set paths
DATADIR=/rds/general/project/lms-ware-raw/live/internal/sequencing/IGFQ001581_ware_4-5-2023_exome_BHFYLKDSX7_2023-07-20/BHFYLKDSX7
FILELISTSDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/wes_ref_controls/sarcoidosis/
mkdir -p $OUTDIR
REFGEN=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta

SAMPLEID=$(sed -n ${PBS_ARRAY_INDEX}p $FILELISTSDIR/samples$cohort.txt)

#concatenate reads
zcat $DATADIR/*/16_NA/$SAMPLEID/${SAMPLEID}_*_R1_* | bgzip -@ 32 > $OUTDIR/$SAMPLEID.pe.1.fastq.gz
zcat $DATADIR/*/16_NA/$SAMPLEID/${SAMPLEID}_*_R2_* | bgzip -@ 32 > $OUTDIR/$SAMPLEID.pe.2.fastq.gz

fastqheader=$(zcat $OUTDIR/$SAMPLEID.pe.1.fastq.gz | head -n 1)
ID=$(echo $fastqheader | cut -d':' -f1,2,3 | cut -c2-)

bwa mem -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID \
-t 38 $REFGEN $OUTDIR/$SAMPLEID.pe.1.fastq.gz $OUTDIR/$SAMPLEID.pe.2.fastq.gz \
| samtools sort --threads 24 -o $OUTDIR/$SAMPLEID.bam

samtools index --threads 16 $OUTDIR/$SAMPLEID.bam

rm $OUTDIR/$SAMPLEID*fastq.gz