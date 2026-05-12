#!/bin/bash

##Job settings
#PBS -V
#PBS -N lpy_lon_UV

##Job configuration
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -lwalltime=1:30:00
#PBS -J 1-1805

##Output log configuration
#PBS -j oe
#PBS -o /rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/logfiles/

#SET COHORT TO ANALYSE
#########################
cohort=HVOL
centre=LA #enter LA (London/Aswan), F (Florence) or P (Prague)
#########################

mkdir -p $EPHEMERAL/$cohort

#activate anaconda module
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

#activate conda env
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/alignment

#set paths
FILELISTSDIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/filelists
DUMPDIR=/rds/general/user/fmazzaro/ephemeral/$cohort
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

    samtools fastq --threads 16 --reference $REFGEN_OLD -1 $DUMPDIR/$SAMPLEID.pe.unsorted.1.fastq -2 $DUMPDIR/$SAMPLEID.pe.unsorted.2.fastq $FPATH
    samtools sort --threads 16 -n $DUMPDIR/$SAMPLEID.pe.unsorted.1.fastq -o $DUMPDIR/$SAMPLEID.pe.1.fastq
    samtools sort --threads 16 -n $DUMPDIR/$SAMPLEID.pe.unsorted.2.fastq -o $DUMPDIR/$SAMPLEID.pe.2.fastq
    rm $DUMPDIR/$SAMPLEID.pe.unsorted.1.fastq
    rm $DUMPDIR/$SAMPLEID.pe.unsorted.2.fastq

    fastqheader=$(head -n 1 $DUMPDIR/$SAMPLEID.pe.1.fastq)
    ID=$(echo $fastqheader | cut -d':' -f1,2,3,4 | cut -c2-)

    echo "Aligning reads of sample $SAMPLEID..."
    echo $(date)

    bwa mem -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID -t 4 $REFGEN_NEW $DUMPDIR/$SAMPLEID.pe.1.fastq $DUMPDIR/$SAMPLEID.pe.2.fastq \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view --threads 2 -S -b - \
    | samtools sort --threads 2 > $DUMPDIR/$SAMPLEID.bam

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
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view --threads 2 -S -b - \
    | samtools sort --threads 2 > $DUMPDIR/$SAMPLEID.bam

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
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view --threads 2 -S -b - \
    | samtools sort --threads 2 > $DUMPDIR/$SAMPLEID.bam
fi



echo "Working on alignment output of sample $SAMPLEID..."
echo $(date)

samtools view --threads 16 -b -F 1294 $DUMPDIR/$SAMPLEID.bam > $DUMPDIR/$SAMPLEID.discordants.unsorted.bam

samtools view --threads 4 -h $DUMPDIR/$SAMPLEID.bam \
| python3 /rds/general/user/fmazzaro/home/miniforge3/envs/lumpy/share/lumpy-sv-0.3.1-3/scripts/extractSplitReads_BwaMem -i stdin \
| samtools view --threads 4 -Sb - > $DUMPDIR/$SAMPLEID.splitters.unsorted.bam

samtools sort --threads 16 $DUMPDIR/$SAMPLEID.discordants.unsorted.bam -o $DUMPDIR/$SAMPLEID.discordants.bam
samtools sort --threads 16 $DUMPDIR/$SAMPLEID.splitters.unsorted.bam -o $DUMPDIR/$SAMPLEID.splitters.bam

echo "Performing lumpy analysis of sample $SAMPLEID..."
echo $(date)

#activate conda env
conda activate lumpy

#run lumpy-sv
lumpyexpress \
-B $DUMPDIR/$SAMPLEID.bam \
-S $DUMPDIR/$SAMPLEID.splitters.bam \
-D $DUMPDIR/$SAMPLEID.discordants.bam \
-o $DUMPDIR/$SAMPLEID.vcf

#remove BND variants
sed -i '/SVTYPE=BND/d' $DUMPDIR/$SAMPLEID.vcf

#index bam for svtyper to run
samtools index -@ 16 $DUMPDIR/$SAMPLEID.bam

#call genotypes
svtyper \
-B $DUMPDIR/$SAMPLEID.bam \
-i $DUMPDIR/$SAMPLEID.vcf \
-o $DUMPDIR/$SAMPLEID.gt.vcf

#remove reference sites
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/bcftools
bcftools filter --include 'AC > 0' $DUMPDIR/$SAMPLEID.gt.vcf > $DUMPDIR/$SAMPLEID.gt.nonref.vcf

#activate env for annotation
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/vep105

#annotate results
vep \
-i $DUMPDIR/$SAMPLEID.gt.nonref.vcf \
--cache \
--dir_cache $vep_cache_dir \
--species homo_sapiens \
--assembly GRCh38 \
--hgvs \
--canonical \
--verbose \
--force_overwrite \
--vcf \
--af_gnomad \
--allele_number \
-o $DUMPDIR/$SAMPLEID.gt.nonref.annot.vcf

#filter to the 23 analysed genes
filter_vep \
--input_file $DUMPDIR/$SAMPLEID.gt.nonref.annot.vcf \
--format vcf \
--force_overwrite \
--filter "((Feature matches ENST00000290378 \
or Feature matches ENST00000366578 \
or Feature matches ENST00000369085 \
or Feature matches ENST00000533783 \
or Feature matches ENST00000373960 \
or Feature matches ENST00000379802 \
or Feature matches ENST00000368300 \
or Feature matches ENST00000545968 \
or Feature matches ENST00000355349 \
or Feature matches ENST00000228841 \
or Feature matches ENST00000334785 \
or Feature matches ENST00000357525 \
or Feature matches ENST00000369519 \
or Feature matches ENST00000333535 \
or Feature matches ENST00000232975 \
or Feature matches ENST00000344887 \
or Feature matches ENST00000367318 \
or Feature matches ENST00000403994 \
or Feature matches ENST00000589042 \
or Feature matches ENST00000211998 \
or Feature matches ENST00000372980 \
or Feature matches ENST00000292327 \
or Feature matches ENST00000374272) \
and (gnomAD_AF < 0.0001 or not gnomAD_AF) \
and (FILTER is PASS))" \
--output_file $DUMPDIR/$SAMPLEID.gt.nonref.annot.genes.vcf

#remove TTN low-psi variants
#1-extract all TTN variants
bgzip -f $DUMPDIR/$SAMPLEID.gt.nonref.annot.genes.vcf
tabix -p vcf $DUMPDIR/$SAMPLEID.gt.nonref.annot.genes.vcf.gz
tabix -h $DUMPDIR/$SAMPLEID.gt.nonref.annot.genes.vcf.gz chr2:178525989-178830802 > $DUMPDIR/$SAMPLEID.ttn.temp.vcf
bgzip -d $DUMPDIR/$SAMPLEID.gt.nonref.annot.genes.vcf.gz

conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/bedtools
#2-extract all non-TTN variants
bedtools intersect -a $DUMPDIR/$SAMPLEID.gt.nonref.annot.genes.vcf -b $DUMPDIR/$SAMPLEID.ttn.temp.vcf -v -header > $DUMPDIR/$SAMPLEID.nottn.temp.vcf
#3-extract TTN-highpsi variants
bedtools intersect -u -header -a $DUMPDIR/$SAMPLEID.ttn.temp.vcf -b $ttn_exons > $DUMPDIR/$SAMPLEID.ttnhighpsi.temp.vcf
#4-concatenate no-TTN and TTN-highpsi variants
conda activate /rds/general/user/fmazzaro/home/anaconda3/envs/bcftools
bgzip -f $DUMPDIR/$SAMPLEID.ttnhighpsi.temp.vcf
bgzip -f $DUMPDIR/$SAMPLEID.nottn.temp.vcf
tabix -p vcf $DUMPDIR/$SAMPLEID.ttnhighpsi.temp.vcf.gz
tabix -p vcf $DUMPDIR/$SAMPLEID.nottn.temp.vcf.gz
bcftools concat --allow-overlaps $DUMPDIR/$SAMPLEID.ttnhighpsi.temp.vcf.gz $DUMPDIR/$SAMPLEID.nottn.temp.vcf.gz > $DUMPDIR/$SAMPLEID.final.unsorted.vcf
bcftools sort --output-file $DUMPDIR/$SAMPLEID.final.vcf.gz --output-type z $DUMPDIR/$SAMPLEID.final.unsorted.vcf
tabix -p vcf $DUMPDIR/$SAMPLEID.final.vcf.gz

#remove files with no variants
if [[ $(zgrep -v "^#" $DUMPDIR/$SAMPLEID.final.vcf.gz | wc -l) -eq 0 ]]; then
    rm $DUMPDIR/$SAMPLEID.final.vcf.gz $DUMPDIR/$SAMPLEID.final.vcf.gz.tbi
fi

rm $DUMPDIR/$SAMPLEID.bam $DUMPDIR/$SAMPLEID.bam.bai $DUMPDIR/$SAMPLEID.vcf $DUMPDIR/$SAMPLEID.discordants.bam $DUMPDIR/$SAMPLEID.splitters.bam $DUMPDIR/$SAMPLEID.pe.1.fastq $DUMPDIR/$SAMPLEID.pe.2.fastq $DUMPDIR/$SAMPLEID.gt.* $DUMPDIR/$SAMPLEID.nottn.* $DUMPDIR/$SAMPLEID.ttn* $DUMPDIR/$SAMPLEID.final.unsorted.vcf $DUMPDIR/$SAMPLEID.discordants.unsorted.bam $DUMPDIR/$SAMPLEID.splitters.unsorted.bam

echo "Analysis of $SAMPLEID finished!!!"
echo $(date)