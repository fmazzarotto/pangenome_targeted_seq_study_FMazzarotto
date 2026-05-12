#!/bin/bash

##Job settings
#PBS -V
#PBS -N p3_HG002

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=6:00:00

#SET COHORT TO ANALYSE
#########################
cohort=HG002    #NOTE: cohort is called SARCWES_forHG002 in part 1 script - OK to have them different
#########################

echo "Part 3 of GATK gCNV analysis of HG002 and 167 control WES (+HG005) samples has started..."
echo $(date)

#set max n of cpus for JAVA/GATK to use, to avoid PBS killing of the job
export GATK_JAVA_OPTS="-XX:ActiveProcessorCount=16" 

CONTAINERDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/docker_images
TEMPOUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKgCNV/tmp
PLOIDYPRIORS=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gatkgcnv/contig_ploidy_prior_hg38.tsv
INTERVALS=$TEMPOUTDIR/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid_padded_GATKgCNVanalysis.interval_list
FILT_INTERVALS=$TEMPOUTDIR/${cohort}_cohort.gc.filtered.interval_list
ANNOT_INTERVALS=$TEMPOUTDIR/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid_padded_GATKgCNVanalysis.annotated_intervals.txt
REFGENOME=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta

#create appropriate -I flag content for the following commands
INPUTSTRING=$(find "$TEMPOUTDIR" -maxdepth 1 -name "*.tsv" -type f -print0 | xargs -0 printf '-I %s ')

#annotate intervals
if [ ! -e $ANNOT_INTERVALS_RAW ]
then
    singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
    gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" AnnotateIntervals \
    -L $INTERVALS \
    -R $REFGENOME \
    -imr OVERLAPPING_ONLY \
    -O $ANNOT_INTERVALS
fi

#filter intervals, re-annotate them, scatter them in 20 parts and create list of the 20 files
if [ ! -e $TEMPOUTDIR/${cohort}_GATKgCNV_bedfiles_scatter.txt ]
then
    singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
    gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" FilterIntervals \
    -L $INTERVALS \
    --annotated-intervals $ANNOT_INTERVALS \
    $INPUTSTRING \
    -imr OVERLAPPING_ONLY \
    -O $FILT_INTERVALS

    grep '^@' $FILT_INTERVALS > $TEMPOUTDIR/header.txt
    grep -v '^@' $FILT_INTERVALS > $TEMPOUTDIR/body.txt
    mkdir -p $TEMPOUTDIR/bed_shards
    split -l 10000 $TEMPOUTDIR/body.txt $TEMPOUTDIR/bed_shards/shard_

    for f in $TEMPOUTDIR/bed_shards/shard_*; do
        cat $TEMPOUTDIR/header.txt "$f" > "${f}.interval_list"
        rm $f
    done

    #singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
    #gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" SplitIntervals \
    #-L $FILT_INTERVALS \
    #-R $REFGENOME \
    #--scatter-count 20 \
    #--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
    #-O $TEMPOUTDIR/${cohort}_cohort.gc.filtered.scattered.interval_list

    #readlink -f $TEMPOUTDIR/${cohort}_cohort.gc.filtered.scattered.interval_list/* > $TEMPOUTDIR/${cohort}_GATKgCNV_bedfiles_scatter.txt
    readlink -f $TEMPOUTDIR/bed_shards/* > $TEMPOUTDIR/${cohort}_GATKgCNV_bedfiles_scatter.txt

fi

#determine samples' ploidy to create the model
singularity exec --bind /rds:/rds $CONTAINERDIR/gatk_latest.sif \
gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" DetermineGermlineContigPloidy \
-L $TEMPOUTDIR/${cohort}_cohort.gc.filtered.interval_list \
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