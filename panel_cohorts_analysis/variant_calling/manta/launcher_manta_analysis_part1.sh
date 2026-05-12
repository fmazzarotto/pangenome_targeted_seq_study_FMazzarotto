#!/bin/bash

##Job settings
#PBS -V
#PBS -N p1_mnt_DCM

##Job configuration
#PBS -l select=1:ncpus=10:mem=16gb
#PBS -lwalltime=1:00:00
#PBS -J 1-928

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=DCM #set this to DCM, HCM or HVOL
#########################

#activate anaconda module and env
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate manta

#set file to analyze based on the array index
BAMFILE=$(sed -n ${PBS_ARRAY_INDEX}p /rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists/files$cohort.txt)
SAMPLEID=$(basename $BAMFILE .bam)

#set other useful folder paths
DATADIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/$cohort
REFGENOME=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
VEP_CACHE_DIR=/rds/general/project/lms-ware-analysis/live/VEP

if [[ "$cohort" == "DCM" || "$cohort" == "HVOL" ]]
then
    TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38_only${cohort}genes_TTNhighPSI.bed.gz
elif [[ "$cohort" == "HCM" ]]
then
    TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38_only${cohort}genes.bed.gz
else
    echo "Error: Check cohort name."
    exit 1
fi

OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/$SAMPLEID
mkdir -p $OUTDIR

#WARNING: launch this without the --target flag otherwise variants spanning the target but with POS outside it are missed
configManta.py \
--bam $BAMFILE \
--exome \
--referenceFasta $REFGENOME \
--runDir $OUTDIR

#for some reason, manta needs me to create these two folders to avoid failing
mkdir -p $OUTDIR/workspace/svLocusGraph.bin.tmpdir
mkdir -p $OUTDIR/workspace/svHyGen

#run analysis
python $OUTDIR/runWorkflow.py \
--jobs 8 \
--memGb 14

#move relevant file and remove output directory
mv $OUTDIR/results/variants/diploidSV.vcf.gz /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_anyFILTER.vcf.gz
rm -rf $OUTDIR

conda deactivate

#index file and keep only PASS variants
conda activate bcftools
tabix -p vcf /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_anyFILTER.vcf.gz

bcftools view \
-f PASS \
--regions-file $TARGET \
/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_anyFILTER.vcf.gz \
-O z \
-o /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_noannot.vcf.gz

rm /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_anyFILTER.vcf.gz*
conda deactivate

#annotate it and index it
conda activate vep105

vep \
-i /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_noannot.vcf.gz \
--format vcf \
--cache \
--dir_cache $VEP_CACHE_DIR \
--species homo_sapiens \
--assembly GRCh38 \
--hgvs \
--canonical \
--verbose \
--force_overwrite \
--vcf \
--compress_output bgzip \
--af_gnomad \
--allele_number \
-o /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_annot.vcf.gz

tabix -p vcf /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_annot.vcf.gz

rm /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_noannot.vcf.gz
rm /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort/${SAMPLEID}_annot.vcf.gz_summary.html