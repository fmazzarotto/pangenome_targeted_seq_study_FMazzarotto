#!/bin/bash

##Job settings
#PBS -V
#PBS -N gatkhc_HG002

##Job configuration
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -lwalltime=4:00:00

echo "GATKhc analysis of HG002 has started..."
echo $(date)

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

#set max n of cpus for JAVA/GATK to use, to avoid PBS killing of the job
export GATK_JAVA_OPTS="-XX:ActiveProcessorCount=62" 

#set paths
HG002_FOLDER=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/GIAB_SRA_WES_samples
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/GATKhc
REFGEN=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/exome-intersect_v5covered_highconfHG002.bed

mkdir -p $OUTDIR

#align reads
conda activate alignment

echo "Aligning reads..."

echo "Alignment started..."
echo $(date)
if [ ! -e $HG002_FOLDER/HG002.bam.bai ]
then
    bwa mem -R @RG\\tPL:ILLUMINA\\tID:SRR2962669\\tSM:HG002 -t 16 $REFGEN $HG002_FOLDER/HG002_SRR2962669_1.fastq.gz $HG002_FOLDER/HG002_SRR2962669_2.fastq.gz | samtools sort --threads 40 -o $HG002_FOLDER/HG002.bam
    samtools index --threads 24 $HG002_FOLDER/HG002.bam
fi
echo "Alignment ended!"
echo $(date)
exit

#activate conda env
conda activate gatk4.1

#hg38 alignment has already been performed (for ED analysis) - just proceed from the MarkDup step
echo "Marking duplicate reads..."

if [ ! -e $HG002_FOLDER/HG002.md.bam ]
then
    gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" MarkDuplicates \
    --INPUT $HG002_FOLDER/HG002.bam \
    --OUTPUT $HG002_FOLDER/HG002.md.bam \
    --METRICS_FILE $HG002_FOLDER/md.metrics.txt
fi

echo "Performing base quality score recalibration..."

if [ ! -e $HG002_FOLDER/bqsr.matrix.txt ]
then
    gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" BaseRecalibrator \
    -I $HG002_FOLDER/HG002.md.bam \
    --known-sites /rds/general/project/lms-ware-raw/live/resources/known_sites_38/dbSNP156.GRCh38.vcf.gz \
    --known-sites /rds/general/project/lms-ware-raw/live/resources/known_sites_38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /rds/general/project/lms-ware-raw/live/resources/known_sites_38/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -O $HG002_FOLDER/bqsr.matrix.txt \
    -R $REFGEN
fi

if [ ! -e $HG002_FOLDER/HG002.md.bqsr.bam ]
then
    gatk --java-options "-Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyBQSR \
    --bqsr-recal-file $HG002_FOLDER/bqsr.matrix.txt \
    -I $HG002_FOLDER/HG002.md.bam \
    -O $HG002_FOLDER/HG002.md.bqsr.bam
fi

echo "Calling variants..."

if [ ! -e $OUTDIR/HG002_gatkhc_raw.vcf ]
then
    gatk --java-options "-Xmx28G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$EPHEMERAL -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" HaplotypeCaller \
    -I $HG002_FOLDER/HG002.md.bqsr.bam \
    -O $OUTDIR/HG002_gatkhc_raw.vcf \
    -R $REFGEN \
    -L $TARGET \
    --pair-hmm-implementation FASTEST_AVAILABLE \
    --native-pair-hmm-threads 16
fi

if [ ! -e $HG002_FOLDER/HG002_gatkhc_INDEL.vcf.gz ]
then
    gatk VariantRecalibrator \
    -R $REFGEN \
    -V $OUTDIR/HG002_gatkhc_raw.vcf \
    --max-gaussians 4 \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 /rds/general/project/lms-ware-raw/live/resources/known_sites_38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /rds/general/project/lms-ware-raw/live/resources/known_sites_38/dbSNP156.GRCh38.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode INDEL \
    -O $HG002_FOLDER/output_INDEL.recal \
    --tranches-file $HG002_FOLDER/output_INDEL.tranches

    gatk ApplyVQSR \
    -R $REFGEN \
    -V $OUTDIR/HG002_gatkhc_raw.vcf \
    -O $HG002_FOLDER/HG002_gatkhc_INDEL.vcf.gz \
    --truth-sensitivity-filter-level 99.9 \
    --tranches-file $HG002_FOLDER/output_INDEL.tranches \
    --recal-file $HG002_FOLDER/output_INDEL.recal \
    --create-output-variant-index true \
    -mode INDEL

fi

if [ ! -e $OUTDIR/HG002_gatkhc_multiall.vcf.gz ]
then
    gatk VariantRecalibrator \
    -R $REFGEN \
    -V $HG002_FOLDER/HG002_gatkhc_INDEL.vcf.gz \
    --max-gaussians 4 \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /rds/general/project/lms-ware-raw/live/resources/known_sites_38/dbSNP156.GRCh38.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /rds/general/project/lms-ware-raw/live/resources/known_sites_38/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 /rds/general/project/lms-ware-raw/live/resources/known_sites_38/1000G_omni2.5.hg38.vcf.gz \
    --resource:hiconf1kg,known=false,training=true,truth=false,prior=10.0 /rds/general/project/lms-ware-raw/live/resources/known_sites_38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -an QD -an MQ -an FS -an SOR -an MQRankSum -an ReadPosRankSum  \
    -mode SNP \
    -O $HG002_FOLDER/output_SNP.recal \
    --tranches-file $HG002_FOLDER/output_SNP.tranches

    gatk ApplyVQSR \
    -R $REFGEN \
    -V $HG002_FOLDER/HG002_gatkhc_INDEL.vcf.gz \
    -O $OUTDIR/HG002_gatkhc_multiall.vcf.gz \
    --truth-sensitivity-filter-level 99.9 \
    --tranches-file $HG002_FOLDER/output_SNP.tranches \
    --recal-file $HG002_FOLDER/output_SNP.recal \
    --create-output-variant-index true \
    -mode SNP

fi

if [ ! -e $OUTDIR/HG002_gatkhc_min20bp_FINAL.vcf.gz.tbi ]
then
    conda activate bcftools
    
    bcftools norm \
    -m -both $OUTDIR/HG002_gatkhc_multiall.vcf.gz \
    -O z -o $OUTDIR/HG002_gatkhc_FINAL.vcf.gz

    bcftools view \
    --threads 8 \
    --include 'strlen(REF) != strlen(ALT) & (strlen(REF) >= 20 | strlen(ALT) >= 20) & FILTER="PASS"' \
    $OUTDIR/HG002_gatkhc_FINAL.vcf.gz \
    | bcftools sort \
    -O z \
    -o $OUTDIR/HG002_gatkhc_min20bp_FINAL.vcf.gz

    bcftools view \
    --threads 8 \
    --include '(strlen(REF) < 20 & strlen(ALT) < 20) & FILTER="PASS"' \
    $OUTDIR/HG002_gatkhc_FINAL.vcf.gz \
    | bcftools sort \
    -O z \
    -o $OUTDIR/HG002_gatkhc_max19bp_FINAL.vcf.gz

    tabix -p vcf $OUTDIR/HG002_gatkhc_min20bp_FINAL.vcf.gz
    tabix -p vcf $OUTDIR/HG002_gatkhc_max19bp_FINAL.vcf.gz

fi



echo "GATKhc analysis of HG002 finished!!!"
echo $(date)

