#!/bin/bash

##Job settings
#PBS -V
#PBS -N manta_HG002

##Job configuration
#PBS -l select=1:ncpus=16:mem=24gb
#PBS -l walltime=06:00:00

#activate anaconda module and env
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate manta

HG002_BAM="/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/GIAB_SRA_WES_samples/HG002.bam"
REFGENOME="/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta"
OUTDIR="/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/manta"

#did not specify a target file here as Manta would call only variants with POS inside it (even is the variants spans it)

configManta.py \
--bam $HG002_BAM \
--exome \
--referenceFasta $REFGENOME \
--runDir $OUTDIR

#for some reason, manta needs me to create these two folders to avoid failing
mkdir -p $OUTDIR/workspace/svLocusGraph.bin.tmpdir
mkdir -p $OUTDIR/workspace/svHyGen

python $OUTDIR/runWorkflow.py \
--jobs 16 \
--memGb 24

conda deactivate

conda activate bcftools

bcftools view \
-f PASS \
$OUTDIR/results/variants/diploidSV.vcf.gz \
| bcftools norm \
-m -both \
--output-type z \
--output $OUTDIR/results/variants/HG002_manta_min15bp_FINAL.vcf.gz

tabix -p vcf $OUTDIR/results/variants/HG002_manta_min15bp_FINAL.vcf.gz

conda deactivate