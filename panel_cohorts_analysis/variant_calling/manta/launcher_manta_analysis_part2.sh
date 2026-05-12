#!/bin/bash

##Job settings
#PBS -V
#PBS -N p2_mnt

##Job configuration
#PBS -l select=1:ncpus=8:mem=8gb
#PBS -lwalltime=1:00:00

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=HVOL
#########################

#activate anaconda module and env
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate bcftools

#set folder paths and increase max n of files that can be open simultaneously
DATADIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/manta_analysis/$cohort
TABLEIZEDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/loftee/src #location of the tableize.py executable from the LOFTEE VEP plugin
ulimit -n 8000

if [[ "$cohort" == "DCM" || "$cohort" == "HVOL" ]]
then
    TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_2bpOverhang_hg38_only${cohort}genes_TTNhighPSI.bed
elif [[ "$cohort" == "HCM" ]]
then
    TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_2bpOverhang_hg38_only${cohort}genes.bed
else
    echo "Error: Check cohort name."
    exit 1
fi


#merge all Manta variant calls (previously filtered for PASS, target, high-PSI TTN exons, protein-altering effect)
bcftools merge \
--missing-to-ref \
--threads 8 \
--output $DATADIR/manta_${cohort}_all.vcf.gz \
--output-type z \
--write-index \
"$DATADIR"/*.vcf.gz

#conda deactivate

conda activate bcftools
#bcftools view -i 'CSQ!~"intron_variant"' $DATADIR/manta_${cohort}_all.vcf.gz -W -Oz -o $DATADIR/manta_${cohort}_nointronic.vcf.gz
bcftools view -e 'SVTYPE="BND"' --regions-file $TARGET $DATADIR/manta_${cohort}_all.vcf.gz -Oz -o $DATADIR/manta_${cohort}_nobnd.vcf.gz
#conda deactivate

#keep only protein-altering and splice site variants
conda activate vep105
filter_vep \
--input_file $DATADIR/manta_${cohort}_nobnd.vcf.gz \
--format vcf \
--force_overwrite \
--filter "(Consequence is not intron_variant) and (Consequence is not 3_prime_utr_variant) and (Consequence is not 5_prime_utr_variant)" \
| bgzip > $DATADIR/manta_${cohort}_FINAL.vcf.gz
#conda deactivate

#tableize final set of variants
conda activate python2.7
$TABLEIZEDIR/tableize_vcf.py \
--vcf $DATADIR/manta_${cohort}_FINAL.vcf.gz \
--info SVLEN,AC,AC_Het,NS \
--vep_info SYMBOL,Gene,Feature,Consequence,HGVSc,HGVSp,gnomAD_AF \
--canonical_only \
--samples \
--output $DATADIR/manta_${cohort}_FINAL.tsv

#clean folder
rm $DATADIR/manta_${cohort}_nobnd.vcf.gz $DATADIR/manta_${cohort}_FINAL.tsv.log 