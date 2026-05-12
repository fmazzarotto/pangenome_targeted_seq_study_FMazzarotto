#!/bin/bash

##Job settings
#PBS -V
#PBS -N p5_HG002

##Job configuration
#PBS -l select=1:ncpus=32:mem=32gb
#PBS -lwalltime=1:00:00

#set max n of cpus for JAVA/GATK to use, to avoid PBS killing of the job
export GATK_JAVA_OPTS="-XX:ActiveProcessorCount=32" 

WORKDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKgCNV/tmp
CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
REFDICT=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.dict
TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid.bed
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKgCNV

#determine which index was assigned by the pipeline to HG002
HG002_INDEX=$(grep -r "HG002" $WORKDIR/variant-calling-shard1-calls/*/sample_name.txt | grep -oP '(?<=SAMPLE_)\d+')

CALLS_STRING=$(for i in $(seq 1 21); do
  echo "--calls-shard-path $WORKDIR/variant-calling-shard${i}-calls"
done | tr '\n' ' ')

MODEL_STRING=$(for i in $(seq 1 21); do
  echo "--model-shard-path $WORKDIR/variant-calling-shard${i}-model"
done | tr '\n' ' ')

singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" PostprocessGermlineCNVCalls \
${CALLS_STRING} \
${MODEL_STRING} \
--allosomal-contig chrX \
--allosomal-contig chrY \
--contig-ploidy-calls $WORKDIR/ploidy-calls \
--sample-index ${HG002_INDEX} \
--output-genotyped-intervals $OUTDIR/HG002.intervals.vcf.gz \
--output-genotyped-segments $OUTDIR/HG002.segments.vcf.gz \
--output-denoised-copy-ratios $OUTDIR/HG002.copyratios.vcf.gz \
--sequence-dictionary $REFDICT

#keep only variant sites in the target region
#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate bcftools

#this excludes variants called as reference sites and no calls
bcftools view \
--regions-file $TARGET \
--include 'GT!="0/0" & GT!="0" & GT!="./." & GT!="."' \
$OUTDIR/HG002.segments.vcf.gz \
| bcftools sort \
--output-type z \
--output-file $OUTDIR/GATKgCNV_HG002_FINAL.vcf.gz

tabix -p vcf $OUTDIR/GATKgCNV_HG002_FINAL.vcf.gz

#clean tmp dir
find $WORKDIR -mindepth 1 -depth -delete
rm -r $WORKDIR
