#!/bin/bash

##Job settings
#PBS -V
#PBS -N lumpy_HG002

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=10:00:00

echo "Lumpy-sv analysis of HG002 has started..."
echo $(date)

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

#activate conda env
conda activate alignment

#set paths
HG002_FOLDER=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/GIAB_SRA_WES_samples
OUTDIR=/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/lumpy
REFGEN=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta
TARGET=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/exome-intersect_v5covered_highconfHG002.bed

mkdir -p $OUTDIR

#realign sample following lumpy-sv best practices
ID=SRR2962669
SAMPLEID=HG002

if [ ! -e  $HG002_FOLDER/HG002_lumpy.bam ]
then
    bwa mem -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID \
    -t 8 $REFGEN $HG002_FOLDER/HG002_SRR2962669_1.fastq.gz $HG002_FOLDER/HG002_SRR2962669_2.fastq.gz \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view --threads 4 -S -b - \
    | samtools sort --threads 4 > $HG002_FOLDER/HG002_lumpy.bam
fi

if [ ! -e  $HG002_FOLDER/HG002_lumpy.bam.bai ]
then
    samtools index $HG002_FOLDER/HG002_lumpy.bam
fi

#create files needed by lumpy
if [ ! -e  $HG002_FOLDER/HG002_discordants.unsorted.bam ]
then
    samtools view --threads 16 -b -F 1294 $HG002_FOLDER/HG002_lumpy.bam > $HG002_FOLDER/HG002_discordants.unsorted.bam
fi

if [ ! -e  $HG002_FOLDER/HG002_splitters.unsorted.bam ]
then
    samtools view --threads 8 -h $HG002_FOLDER/HG002_lumpy.bam \
    | python3 /rds/general/user/fmazzaro/home/miniforge3/envs/lumpy/share/lumpy-sv-0.3.1-3/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view --threads 8 -Sb - > $HG002_FOLDER/HG002_splitters.unsorted.bam
fi

if [ ! -e  $HG002_FOLDER/HG002_discordants.bam ]
then
    samtools sort --threads 16 $HG002_FOLDER/HG002_discordants.unsorted.bam -o $HG002_FOLDER/HG002_discordants.bam
fi

if [ ! -e  $HG002_FOLDER/HG002_splitters.bam ]
then
    samtools sort --threads 16 $HG002_FOLDER/HG002_splitters.unsorted.bam -o $HG002_FOLDER/HG002_splitters.bam
fi

rm $HG002_FOLDER/*unsorted*

#activate conda env
conda activate lumpy

#run lumpy-sv
if [ ! -e  $OUTDIR/HG002_lumpy_raw.vcf ]
then
    lumpyexpress \
    -B $HG002_FOLDER/HG002_lumpy.bam \
    -S $HG002_FOLDER/HG002_splitters.bam \
    -D $HG002_FOLDER/HG002_discordants.bam \
    -o $OUTDIR/HG002_lumpy_raw.vcf
fi

if [ ! -e  $OUTDIR/HG002_lumpy_raw_fixed_sorted.vcf ]
then
    #manually add contigs to the header (to enable sorting and following steps)
    awk '{print "##contig=<ID="$1",length="$2">"}' $REFGEN.fai > $OUTDIR/contig_lines.txt
    grep -E '^##' $OUTDIR/HG002_lumpy_raw.vcf > $OUTDIR/metadata_only.txt
    grep -E '^#CHROM' $OUTDIR/HG002_lumpy_raw.vcf > $OUTDIR/chrom_line.txt
    grep -v '^#' $OUTDIR/HG002_lumpy_raw.vcf > $OUTDIR/data_lines.txt

    cat $OUTDIR/metadata_only.txt $OUTDIR/contig_lines.txt $OUTDIR/chrom_line.txt $OUTDIR/data_lines.txt | bgzip -c > $OUTDIR/HG002_lumpy_raw_fixed.vcf.gz

    bcftools sort $OUTDIR/HG002_lumpy_raw_fixed.vcf.gz -O v -o $OUTDIR/HG002_lumpy_raw_fixed_sorted.vcf

    rm $OUTDIR/*.txt
fi

#call genotypes and exclude BND variants
if [ ! -e  $OUTDIR/HG002_lumpy_gt.vcf ]
then
    svtyper \
    -B $HG002_FOLDER/HG002_lumpy.bam \
    -i $OUTDIR/HG002_lumpy_raw_fixed_sorted.vcf \
    > $OUTDIR/HG002_lumpy_gt.vcf
fi

#filter contigs, remove BND variants
if [ ! -e  $OUTDIR/HG002_lumpy_gt_clean.vcf.gz ]
then
    bcftools view \
    $OUTDIR/HG002_lumpy_gt.vcf \
    -t 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY' \
    -e 'INFO/SVTYPE="BND"' \
    -O z \
    -o $OUTDIR/HG002_lumpy_gt_clean.vcf.gz
fi

#index file
if [ ! -e  $OUTDIR/HG002_lumpy_gt_clean.vcf.gz.tbi ]
then
    tabix -p vcf $OUTDIR/HG002_lumpy_gt_clean.vcf.gz
fi

#remove reference sites and restrict to target
if [ ! -e  $OUTDIR/HG002_lumpy_FINAL.vcf ]
then
    bcftools filter \
    --include 'AC > 0' \
    $OUTDIR/HG002_lumpy_gt_clean.vcf.gz \
    -O v \
    | bcftools norm -m -both | bcftools sort -O z -o $OUTDIR/HG002_lumpy_FINAL.vcf.gz

    tabix -p vcf $OUTDIR/HG002_lumpy_FINAL.vcf.gz
fi

echo "Analysis of $SAMPLEID finished!!!"
echo $(date)