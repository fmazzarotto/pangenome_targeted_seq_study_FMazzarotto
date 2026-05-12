## Files
This repo contains the main scripts used in developing the study described in _"Pangenomes aid accurate detection of large insertion and deletions from targeted sequencing: the case of cardiomyopathies"_ by F Mazzarotto et al.

Folder _"HG002_analysis"_ contains all scripts applied on HG002 raw sequencing data, callers' benchmarking results and the target BED files used in the analysis:
- _alignment_and_variant_calling_: this folder contains all the scripts used to align HG002 raw sequencing reads to the appropriate reference genome (GRCh38-based "global" pangenome reference for Velsera GRAF, GRCh38 reference for other tools), and to perform variant calling. Scripts are named according to the tool used, and include alignment steps. GATK-gCNV scripts are split into 5 parts (HG002_gcnv_part1..part5) to accommodate the SLURM array-based and single-threaded execution structure of the pipeline.
- _callers_benchmarking_: this folder includes the full benchmarking results obtained with Truvari (for variants with a minimum size of 20bp) and RTG Tools vcfeval (for variants in the 1-19bp range) for Manta, LUMPY, Velsera GRAF and GATK HaplotypeCaller, with the 1-19bp comparison performed on the last two. In addition, this folder contains variant calls made by ExomeDepth and GATK-gCNV, which were compared with HG002 truth set manually, using bedtools intersect.
- _HG002_truth_set_: this folder contains a small text file specifying the source of the HG002 VCF file with the full set of ground truth variants, besides the BED file listing high-confidence regions for HG002.
- _reference_cohort_creation_for_ExomeDepth_and_GATKgCNV_: this contains the script used to align exome sequencing reads from 168 samples (GIAB HG005 and 167 internally available individuals) against GRCh38, as this is the sample set used as reference cohort for ExomeDepth and GATK-gCNV.

Folder _"panel_cohorts_analysis"_ contains all scripts used to analyse the 3774 cohort samples featured in the study with Manta, ExomeDepth, GATK-gCNV and LUMPY:
- _alignment_: this folder contains the script used to align raw sequencing data from all samples against GRCh38.
- _variant_calling_: this contains four caller-specific subfolders with executables used to call variants on all samples. Whenever split into multiple scripts, these are numbered sequentially.
- _post_processing_results_: this folder contains scripts used to post-process variant calling results (e.g. removing variants altering low-PSI exons in TTN) and/or to reformat and harmonize outputs for downstream analysis.

Commands used in processing cohort samples with GRAF are the same used on HG002 (and listed in the manuscript's Supplementary Notes). Commands used in processing cohort samples with GATK-HaplotypeCaller are those suggested by GATK Best Practices, using default settings. Software required to reproduce the results is available as Docker images, as indicated in the manuscript's Supplementary Notes.

## Citation
If you use this code in academic work, please cite our manuscript.
