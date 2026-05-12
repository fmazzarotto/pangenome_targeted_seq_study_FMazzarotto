#!/bin/bash

set -e

#activate conda
module load anaconda3/personal
eval "$(conda shell.bash hook)"

#specify paths and relevant files
resources_dir="/rds/general/user/fmazzaro/home/WORK/Large_variants/resources" #folder with resources e.g. TTN exons, executables etc
results_dir="/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/hard_filtered_postqc/hvol" #folder where results are to be written
data_dir="/rds/general/user/fmazzaro/home/WORK/Large_variants/data/hard_filtered_postqc/hvol" #raw data (VCF files)
ttn_exons=$resources_dir/ttnexons/TTN_highPSI_GRCh38.bed #bed file with coordinates of high-PSI TTN exons
tableize_dir=$resources_dir/loftee/src #location of the tableize.py executable from the LOFTEE VEP plugin
vep_cache_dir="/rds/general/project/lms-ware-analysis/live/VEP"


###DESCRIPTION: this script works by filtering/processing the VEP-annotated 'raw' VCF files provided by the Sevenbridges team.

###The --filter field of the filter_vep command at the beginning serves to filter by gene(s), so it needs to be modified according to the cohort under analysis.
###The three filter strings to be used are the following:

##############BY TRANSCRIPT:
###DCM: Feature matches ENST00000290378 or Feature matches ENST00000366578 or Feature matches ENST00000369085 or Feature matches ENST00000373960 or Feature matches ENST00000379802 or Feature matches ENST00000368300 or Feature matches ENST00000355349 or Feature matches ENST00000334785 or Feature matches ENST00000357525 or Feature matches ENST00000369519 or Feature matches ENST00000333535 or Feature matches ENST00000232975 or Feature matches ENST00000344887 or Feature matches ENST00000367318 or Feature matches ENST00000403994 or Feature matches ENST00000589042 or Feature matches ENST00000211998 or Feature matches ENST00000372980
###HCM: Feature matches ENST00000290378 or Feature matches ENST00000366578 or Feature matches ENST00000533783 or Feature matches ENST00000545968 or Feature matches ENST00000355349 or Feature matches ENST00000228841 or Feature matches ENST00000357525 or Feature matches ENST00000232975 or Feature matches ENST00000344887 or Feature matches ENST00000367318 or Feature matches ENST00000403994 or Feature matches ENST00000372980 or Feature matches ENST00000292327 or Feature matches ENST00000374272
###HVOL (DCM+HCM): Feature matches ENST00000290378 or Feature matches ENST00000366578 or Feature matches ENST00000369085 or Feature matches ENST00000533783 or Feature matches ENST00000373960 or Feature matches ENST00000379802 or Feature matches ENST00000368300 or Feature matches ENST00000545968 or Feature matches ENST00000355349 or Feature matches ENST00000228841 or Feature matches ENST00000334785 or Feature matches ENST00000357525 or Feature matches ENST00000369519 or Feature matches ENST00000333535 or Feature matches ENST00000232975 or Feature matches ENST00000344887 or Feature matches ENST00000367318 or Feature matches ENST00000403994 or Feature matches ENST00000589042 or Feature matches ENST00000211998 or Feature matches ENST00000372980 or Feature matches ENST00000292327 or Feature matches ENST00000374272
######################################################################################################################################


###NOTE THAT BEFORE WRITING THIS SCRIPT I CHANGED RAW VCF FILE NAMES AS FOLLOWS (format: caller_cohort_genome_raw.vcf):
#bwa_gatk_dcm_merged_vep.vep.vcf -> gatk_dcm_l_raw.vcf
#bwa_manta_dcm_vep.vep.vcf -> manta_dcm_l_raw.vcf
#graf_dcm_all_wdcmsvgV2_merged.filtered.norm.vcf -> graf_dcm_d_raw.vcf
#bwa_gatk_hcm_vep.vep.vcf -> gatk_hcm_l_raw.vcf
#bwa_manta_hcm_vep.vep.vcf -> manta_hcm_l_raw.vcf
#graf_hcm_all_whcmsvg_merged.filtered.norm.vcf -> graf_hcm_h_raw.vcf
#bwa_gatk_healthy_vep.vep.vcf -> gatk_hvol_l_raw.vcf
#bwa_manta_healthy_vep.vep.vcf -> manta_hvol_l_raw.vcf
#graf_hvol_all_wdcmsvgV2_merged.filtered.norm.vcf -> graf_hvol_d_raw.vcf
#graf_hvol_all_whcmsvg_merged.filtered.norm.vcf -> graf_hvol_h_raw.vcf

###SAME FOR HARD-FILTERED GRAF VCFs AFTER THE ROUND OF VALIDATIONS
#dcm_dcmsvg.norm.filtered.1_928.norm.merged.vcf -> graf_dcm_d_hf.vcf
#hcm_hcmsvg.norm.filtered.1_1041.norm.merged.vcf -> graf_hcm_h_hf.vcf
#hvol_dcmsvg.norm.filtered.1_1805.norm.merged.vcf -> graf_hvol_d_hf.vcf
#hvol_hcmsvg.norm.filtered.1_1805.norm.merged.vcf -> graf_hvol_h_hf.vcf

###ALL INPUT FILES WERE THEN COMPRESSED WITH bgzip -f


for fullpath_input_file in $data_dir/*.gz
do
	input_file=$(echo $fullpath_input_file | rev | cut -d'/' -f1 | rev)

	prefix=$(cut -d'_' -f1,2,3 <<< "$input_file")


	#reannotate raw files and extract variants with no CSQ field
	if [ ! -e $results_dir/${prefix}_reannotated.vcf ]
	then

		#if files are produced by Manta, remove BND variants (as VEP is unable to parse them and outputs an error [can't detect input format])
		if [[  $input_file == "manta"* ]]
		then
			noncomp_input=$(echo "$input_file" | rev | cut -d. -f2- | rev)
			printf "\n\n\nRemoving BND variants from $input_file...\n\n\n"

			source activate python2.7 #activate conda environment with python 2.7 and tabix (needed for tableize)
			sed '/SVTYPE=BND/d' $data_dir/$noncomp_input > $data_dir/${prefix}_nobnd.vcf #this step filters out 78 (DCM) - 95% (HVOL) variants
			bgzip -f $data_dir/${prefix}_nobnd.vcf
			input_file=${prefix}_nobnd.vcf.gz
		fi

		source activate vep105 #activate conda environment with ensembl-vep=105

		printf "\n\n\n1) Reannotating $input_file with the --allele-number option to allow processing with tableize later on...\n\n\n"

		#reannotate raw files with the allele_number vep option to allow processing with tableize

		vep \
        	-i $data_dir/$input_file \
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
        	-o $results_dir/${prefix}_reannotated.vcf

		printf "\n\n\n2) Extracting no-CSQ variants (ALT=*) from $results_dir/${prefix}_reannotated.vcf...\n\n\n"

		#extract ALT=* variants from the reannotated files (need to do it now, as they would be filtered out by the following commands based on ENST IDs)

		source activate bcftools #activate conda environment with bcftools (v1.9)

		bcftools filter -i 'ALT="*"' $results_dir/${prefix}_reannotated.vcf > $results_dir/${prefix}_noCSQ_raw.vcf

		source activate python2.7 #activate conda environment with python 2.7 and tabix (needed for tableize)

		$tableize_dir/tableize_vcf.py \
		--vcf $results_dir/${prefix}_noCSQ_raw.vcf \
		--info AC,AC_Het \
		--vep_info gnomAD_AF \
		--samples \
		--output $results_dir/${prefix}_noCSQ.tsv

	fi

	#extract rare variants in the relevant genes by filtering by canonical transcript
	if [ ! -e $results_dir/${prefix}_genes.vcf.gz ]
	then
		printf "\n\n\n3) Extracting variants in relevant genes from ${prefix}_reannotated.vcf...\n\n\n"

		source activate vep105 #activate conda environment with ensembl-vep=105

		filter_vep \
		--input_file $results_dir/${prefix}_reannotated.vcf \
		--format vcf \
		--force_overwrite \
		--filter "((Feature matches ENST00000290378 or Feature matches ENST00000366578 or Feature matches ENST00000369085 or Feature matches ENST00000533783 or Feature matches ENST00000373960 or Feature matches ENST00000379802 or Feature matches ENST00000368300 or Feature matches ENST00000545968 or Feature matches ENST00000355349 or Feature matches ENST00000228841 or Feature matches ENST00000334785 or Feature matches ENST00000357525 or Feature matches ENST00000369519 or Feature matches ENST00000333535 or Feature matches ENST00000232975 or Feature matches ENST00000344887 or Feature matches ENST00000367318 or Feature matches ENST00000403994 or Feature matches ENST00000589042 or Feature matches ENST00000211998 or Feature matches ENST00000372980 or Feature matches ENST00000292327 or Feature matches ENST00000374272) \
		and (gnomAD_AF < 0.0001 or not gnomAD_AF) \
		and (FILTER is PASS))" \
		--output_file $results_dir/${prefix}_genes.vcf

		cp $results_dir/${prefix}_genes.vcf $results_dir/${prefix}_genes_anypsi.vcf #this is just needed for both the if -e clause and the for loop below to work
	fi
done

# exclude variants altering only low-PSI TTN exons
for anypsi in $results_dir/*genes_anypsi.vcf
do

	filename=$(echo $anypsi | rev | cut -d'/' -f1 | rev)
	prefix=$(cut -d'_' -f1,2,3 <<< "$filename")

	#remove these files ($anypsi is a copy) as they are recreated in this for loop with only high-PSI exon variants in TTN
	rm -f $results_dir/${prefix}_genes.vcf

	#separate TTN from non-TTN variants
	if [ ! -e $results_dir/${prefix}_nottn.vcf ]
	then
		printf "\n\n\n4) Excluding variants altering only low-PSI TTN exons from $filename...\n\n\n"

		source activate python2.7 #activate conda environment with python 2.7 and tabix (needed for tableize)

		cd $results_dir
		bgzip -f $anypsi
		tabix -p vcf $anypsi.gz
		tabix -h -p vcf $anypsi.gz chr2:178525989-178830802 > ${prefix}_ttn.vcf
		bgzip -d $anypsi.gz

		source activate bedtools #activate conda environment with bedtools
		bedtools intersect \
		-a $anypsi \
		-b ${prefix}_ttn.vcf \
		-v -header > ${prefix}_nottn.vcf

		bedtools intersect \
		-a ${prefix}_ttn.vcf \
		-b $ttn_exons \
		-u -header > ${prefix}_ttnhighpsi.vcf

		source activate python2.7 #activate conda environment with python 2.7 and tabix (needed for tableize)
		bgzip -f ${prefix}_nottn.vcf
		bgzip -f ${prefix}_ttnhighpsi.vcf
		tabix -p vcf ${prefix}_nottn.vcf.gz
		tabix -p vcf ${prefix}_ttnhighpsi.vcf.gz

		source activate bcftools #activate conda environment with bcftools (v1.9)
		bcftools concat \
		--allow-overlaps \
		${prefix}_nottn.vcf.gz \
		${prefix}_ttnhighpsi.vcf.gz > ${prefix}_genes_unsorted.vcf.gz

		bcftools sort \
		--output-file ${prefix}_genes.vcf.gz \
		--output-type z \
		${prefix}_genes_unsorted.vcf.gz

	fi
done


#save short variants and indels in separate files to check later if large variant carriers have short P/LP variants
#commands were established after ensuring that the only variant types in files were snps,indels and mnps (the latter only in graf files and max 3-bases long)
if [ ! -e $results_dir/graf_hcm_h_genes_indels.vcf ]
then
	printf "\n\n\n5) Splitting short and long variants...\n\n\n"

	source activate bcftools #activate conda environment with bcftools (v1.9)
	for input_file in $results_dir/*genes.vcf.gz
	do
		filename=$(echo $input_file | rev | cut -d'/' -f1 | rev)
		prefix=$(cut -d'_' -f1,2,3 <<< "$filename")

		[[ $filename == manta* ]] && continue		#ignore manta files (all variants need to be kept)

		cd $results_dir
		bcftools view \
		--output-type v \
		--output-file ${prefix}_genes_short_variants.vcf \
		--types snps,mnps \
		$input_file

		bcftools view \
		--output-type v \
		--output-file ${prefix}_genes_indels.vcf \
		--types indels \
		$input_file
	done
fi

#extract indels that are min 20 bases long
if [ ! -e $results_dir/graf_hcm_h_FINAL.vcf ]
then
	printf "\n\n\n6) Extracting variants min 20 bases long...\n\n\n"

	source activate gatk #activate conda environment with gatk4

	for input_file in $results_dir/*indels.vcf
	do
		filename=$(echo $input_file | rev | cut -d'/' -f1 | rev)
		prefix=$(cut -d'_' -f1,2,3 <<< "$filename")

		[[ $filename == manta* ]] && continue		#ignore manta files (all variants need to be kept)

		cd $results_dir
		gatk SelectVariants \
		--min-indel-size 20 \
		--variant $input_file \
		--output ${prefix}_FINAL.vcf
	done
fi


###NOTE: commented out lines from here are excluded when processing hard-filtered, latest-version GRAF files (otherwise it throws an error because there are no Manta files)
# #rename manta files
# if [ ! -e $results_dir/manta_hcm_l_FINAL.vcf ]
# then
# 	source activate python2.7 #activate conda environment with python 2.7 and tabix (needed for tableize)
# 	cd $results_dir
# 	bgzip -d -c manta_dcm_l_genes.vcf.gz > manta_dcm_l_FINAL.vcf
# 	bgzip -d -c manta_hcm_l_genes.vcf.gz > manta_hcm_l_FINAL.vcf
# 	bgzip -d -c manta_hvol_l_genes.vcf.gz > manta_hvol_l_FINAL.vcf
# fi

#tableize the final VCFs to convert it into tables
source activate python2.7 #activate conda environment with python 2.7 and tabix (needed for tableize)
if [ ! -e $results_dir/manta_hcm_l_RESULTS.tsv ]
then
	for input_file in $results_dir/*FINAL.vcf
	do

		filename=$(echo $input_file | rev | cut -d'/' -f1 | rev)
		prefix=$(cut -d'_' -f1,2,3 <<< "$filename")

		if [[ "$filename" == *"manta"* ]] #condition needed just for extracting also the SVLEN field from Manta files
		then
			$tableize_dir/tableize_vcf.py \
			--vcf $input_file \
			--info SVLEN,AC,AC_Het,NS \
			--vep_info SYMBOL,Gene,Feature,Consequence,HGVSc,HGVSp,gnomAD_AF \
			--canonical_only \
			--samples \
			--output $results_dir/${prefix}_RESULTS.tsv
		else
			$tableize_dir/tableize_vcf.py \
			--vcf $input_file \
			--info AC,AC_Het,NS \
			--vep_info SYMBOL,Gene,Feature,Consequence,HGVSc,HGVSp,gnomAD_AF \
			--canonical_only \
			--samples \
			--output $results_dir/${prefix}_RESULTS.tsv
		fi
	done
fi

#clean folders
rm -f $results_dir/*.log $results_dir/*.tbi $results_dir/*summary* $results_dir/*warning* $results_dir/*_ttn.* $results_dir/*.idx $results_dir/*ttnhighpsi* $results_dir/*unsorted*
rm -f $data_dir/*nobnd* #remove the nobnd files created in the data_dir

conda deactivate
