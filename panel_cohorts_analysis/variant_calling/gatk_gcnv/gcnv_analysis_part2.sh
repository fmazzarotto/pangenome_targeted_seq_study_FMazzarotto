#!/bin/bash

##Job settings
#PBS -V
#PBS -N gcnv_p2_HVOL

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=36:00:00

#SET COHORT TO ANALYSE
#########################
cohort=DCM
#########################

CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
TEMPOUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/${cohort}_temp
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/${cohort}
mkdir -p $OUTDIR
PLOIDYPRIORS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gatkgcnv/contig_ploidy_prior_hg38.tsv
INTERVALS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.interval_list
ANNOT_INTERVALS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.annotated.tsv
REFDICT=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.dict

#create appropriate -I flag content for the following commands
INPUTSTRING=$(find "$TEMPOUTDIR" -maxdepth 1 -name "*.tsv" -type f -print0 | xargs -0 printf '-I %s ')

singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" FilterIntervals \
-L $INTERVALS \
--annotated-intervals $ANNOT_INTERVALS \
$INPUTSTRING \
-imr OVERLAPPING_ONLY \
-O $OUTDIR/${cohort}_cohort.gc.filtered.interval_list

singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" DetermineGermlineContigPloidy \
-L $OUTDIR/${cohort}_cohort.gc.filtered.interval_list \
-XL chrY:10001-2781479 \
-XL chrX:10001-2781479 \
-XL chrY:56887903-57217415 \
-XL chrX:155701383-156030895 \
--interval-merging-rule OVERLAPPING_ONLY \
$INPUTSTRING \
--contig-ploidy-priors $PLOIDYPRIORS \
--output $TEMPOUTDIR \
--output-prefix ploidy \
--verbosity DEBUG

singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GermlineCNVCaller \
--run-mode COHORT \
-L $OUTDIR/${cohort}_cohort.gc.filtered.interval_list \
$INPUTSTRING \
--contig-ploidy-calls $TEMPOUTDIR/ploidy-calls \
--annotated-intervals $ANNOT_INTERVALS \
--interval-merging-rule OVERLAPPING_ONLY \
--output $TEMPOUTDIR \
--output-prefix variant-calling \
--verbosity DEBUG