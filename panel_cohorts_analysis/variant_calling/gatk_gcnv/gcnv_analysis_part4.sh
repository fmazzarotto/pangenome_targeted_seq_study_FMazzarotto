#!/bin/bash

##Job settings
#PBS -V
#PBS -N p4_HVOL

##Job configuration
#PBS -l select=1:ncpus=4:mem=4gb
#PBS -lwalltime=2:00:00

#activate conda and the R conda env
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

for cohort in DCM HCM
do
    conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/python2.7
    DATADIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/$cohort
    TABLEIZEDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/loftee/src
    SCRIPTSDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts

    for f in $DATADIR/*.vcf.gz
    do

        SAMPLEID=$(echo "$f" | sed 's/.FINAL.vcf.gz$//')

        #mock SVLEN field extraction - will be NA but need to specify at least 1
        $TABLEIZEDIR/tableize_vcf.py \
        --vcf $f \
        --info SVLEN \
        --include_id \
        --samples \
        --output $SAMPLEID.tsv

        rm $SAMPLEID.tsv.log
        
    done

    conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/r431

    for f in $DATADIR/*.tsv
    do

        echo $f
        Rscript $SCRIPTSDIR/filter_ttnexons_gcnv_results.R "$f" "$cohort"

    done

    awk 'NR==1 || (FNR>1 && $1!="NA")' $DATADIR/*_TTNhp.tsv > $DATADIR/FINAL_${cohort}_VARIANTS.tsv

done

conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/r431
Rscript $SCRIPTSDIR/postprocess_gcnv_calls.R