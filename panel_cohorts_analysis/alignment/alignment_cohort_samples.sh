#!/bin/bash

##Job settings
#PBS -V
#PBS -N aln_F_HCM

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=5:00:00
#PBS -J 1-175

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=HCM #enter HCM, DCM or HVOL
centre=F #enter LA (London/Aswan), F (Florence) or P (Prague)
#########################


OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/$cohort
mkdir -p $OUTDIR

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

#activate conda env
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/alignment

#set paths
FILELISTSDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists
REFGEN_OLD=/rds/general/project/lms-ware-raw/live/resources/reference/UCSC_hg19/allchrom.Chr1ToChrM.validated.fa
REFGEN_NEW=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
vep_cache_dir=/rds/general/project/lms-ware-analysis/live/VEP
ttn_exons=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/ttnexons/TTN_highPSI_GRCh38.bed #bed file with coordinates of high-PSI TTN exons

if [ "$centre" = "LA" ]; then
    FPATH=$(sed -n ${PBS_ARRAY_INDEX}p $FILELISTSDIR/files${cohort}_London_Aswan.txt)
    FNAME=$(basename $FPATH)
    SAMPLEID=$(echo $FNAME | cut -d . -f 1)

    echo "Extracting reads of sample $SAMPLEID..."
    echo $(date)

    samtools fastq --threads 16 --reference $REFGEN_OLD -1 $OUTDIR/$SAMPLEID.pe.unsorted.1.fastq -2 $OUTDIR/$SAMPLEID.pe.unsorted.2.fastq $FPATH
    samtools sort --threads 16 -n $OUTDIR/$SAMPLEID.pe.unsorted.1.fastq -o $OUTDIR/$SAMPLEID.pe.1.fastq
    samtools sort --threads 16 -n $OUTDIR/$SAMPLEID.pe.unsorted.2.fastq -o $OUTDIR/$SAMPLEID.pe.2.fastq
    rm $OUTDIR/$SAMPLEID.pe.unsorted.1.fastq
    rm $OUTDIR/$SAMPLEID.pe.unsorted.2.fastq

    fastqheader=$(head -n 1 $OUTDIR/$SAMPLEID.pe.1.fastq)
    ID=$(echo $fastqheader | cut -d':' -f1,2,3,4 | cut -c2-)

    echo "Aligning reads of sample $SAMPLEID..."
    echo $(date)

    bwa mem -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID -t 4 $REFGEN_NEW $EPHEMERAL/$SAMPLEID.pe.1.fastq $EPHEMERAL/$SAMPLEID.pe.2.fastq \
    | samtools view --threads 2 -S -b - \
    | samtools sort --threads 2 > $OUTDIR/$SAMPLEID.bam

    samtools index $OUTDIR/$SAMPLEID.bam

elif [ "$centre" = "F" ]; then
    sed 's/_[^/]*$/_/' $FILELISTSDIR/files${cohort}_Florence.txt | sort -u > $EPHEMERAL/files${cohort}_Florence_prefix_temp.txt
    FPATH=$(sed -n ${PBS_ARRAY_INDEX}p $EPHEMERAL/files${cohort}_Florence_prefix_temp.txt)
    SAMPLEID=$(basename $FPATH _)
    #create R1 and R2 fastq files
    zcat ${FPATH}*_R1_*.fastq.gz | bgzip --threads 16 > $EPHEMERAL/${SAMPLEID}.pe.1.unsorted.fastq.gz
    zcat ${FPATH}*_R2_*.fastq.gz | bgzip --threads 16 > $EPHEMERAL/${SAMPLEID}.pe.2.unsorted.fastq.gz
    samtools sort --threads 16 -n $EPHEMERAL/${SAMPLEID}.pe.1.unsorted.fastq.gz -o $EPHEMERAL/${SAMPLEID}.pe.1.fastq.gz
    samtools sort --threads 16 -n $EPHEMERAL/${SAMPLEID}.pe.2.unsorted.fastq.gz -o $EPHEMERAL/${SAMPLEID}.pe.2.fastq.gz

    fastqheader=$(zless $EPHEMERAL/${SAMPLEID}.pe.1.fastq.gz | head -n1)
    ID=$(echo $fastqheader | cut -d':' -f1,2,3,4 | cut -c2-)

    echo "Aligning reads of sample $SAMPLEID..."
    echo $(date)

    bwa mem -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID -t 4 $REFGEN_NEW $EPHEMERAL/$SAMPLEID.pe.1.fastq.gz $EPHEMERAL/$SAMPLEID.pe.2.fastq.gz \
    | samtools view --threads 2 -S -b - \
    | samtools sort --threads 2 > $OUTDIR/$SAMPLEID.bam

    samtools index $OUTDIR/$SAMPLEID.bam

elif [ "$centre" = "P" ]; then
    sed 's/[_-][^/]*$/_/' $FILELISTSDIR/files${cohort}_Prague.txt | sort -u > $EPHEMERAL/files${cohort}_Prague_prefix_temp.txt
    FPATH=$(sed -n ${PBS_ARRAY_INDEX}p $EPHEMERAL/files${cohort}_Prague_prefix_temp.txt)
    SAMPLEID=$(basename $FPATH _)

    #map variables to the R1 and R2 file for the sample in question
    r1file=$(ls /rds/general/project/lms-ware-raw/live/external/Prague_HCM/data/$SAMPLEID*_R1*.fastq.gz)
    r2file=$(ls /rds/general/project/lms-ware-raw/live/external/Prague_HCM/data/$SAMPLEID*_R2*.fastq.gz)

    samtools sort --threads 16 -n $r1file -o $EPHEMERAL/${SAMPLEID}.pe.1.fastq.gz
    samtools sort --threads 16 -n $r2file -o $EPHEMERAL/${SAMPLEID}.pe.2.fastq.gz

    fastqheader=$(zless $r1file | head -n1)
    ID=$(echo $fastqheader | cut -d':' -f1,2,3,4 | cut -c2-)

    echo "Aligning reads of sample $SAMPLEID..."
    echo $(date)

    bwa mem -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID -t 4 $REFGEN_NEW $EPHEMERAL/${SAMPLEID}.pe.1.fastq.gz $EPHEMERAL/${SAMPLEID}.pe.2.fastq.gz \
    | samtools view --threads 2 -S -b - \
    | samtools sort --threads 2 > $OUTDIR/$SAMPLEID.bam

    samtools index $OUTDIR/$SAMPLEID.bam
fi