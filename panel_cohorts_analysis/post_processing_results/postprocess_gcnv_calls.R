
library(dplyr)

####HVOL####

#import HVOL vars and remove duplicate entries
hvol_vars = read.table("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/HVOL/FINAL_HVOL_VARIANTS.tsv",header=T,sep="\t")
hvol_vars = unique(hvol_vars)

#sort by QUAL
hvol_vars = hvol_vars[order(hvol_vars$QUAL, decreasing=TRUE),]

#using dplyr, keep the top X variants so to have 18 carriers (keep the top X variants so to have 18 carriers)
# Find the minimum number of rows needed to get 18 unique carriers
unique_samples = character(0)
rows_to_keep = 0

for (i in 1:nrow(hvol_vars)) {
  unique_samples <- unique(c(unique_samples, hvol_vars$SAMPLES[i]))
  if (length(unique_samples) >= 9) {
    rows_to_keep <- i
    break
  }
}

hvol_vars_prioritized = hvol_vars[1:rows_to_keep, ]

#derive minimum QUAL to keep a variant
minqual = min(hvol_vars_prioritized$QUAL)

#fix order of columns etc
hvol_vars_prioritized$CHROM = gsub("chr","",hvol_vars_prioritized$CHROM)
hvol_vars_prioritized = cbind(hvol_vars_prioritized[,c(1:6,8:9,7,10:12)],rep("",nrow(hvol_vars_prioritized)),rep("",nrow(hvol_vars_prioritized)),rep("",nrow(hvol_vars_prioritized)),hvol_vars_prioritized[,c(13:21)],rep("HVOL",nrow(hvol_vars_prioritized)),rep("GATKgcnv",nrow(hvol_vars_prioritized)))
colnames(hvol_vars_prioritized)[c(9,13:15,25,26)] = c("GATKgcnv_QUAL","AC","AC_Het","N_SAMPLES","COHORT","CALLER")

#always using dplyr, collapse rows referring to the same variant in different samples
collapsed_hvol_vars = hvol_vars_prioritized %>%
  group_by(across(-c(GATKgcnv_QUAL,SAMPLES))) %>%
  summarise(
    across(everything(), ~paste(., collapse = ",")),
    .groups = "drop"
  )
collapsed_hvol_vars = collapsed_hvol_vars[,colnames(hvol_vars_prioritized)]

write.table(hvol_vars_prioritized,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/HVOL/FINAL_HVOL_VARIANTS_FILTERED.tsv",row.names=F,col.names=T,sep="\t",quote=F)
write.table(collapsed_hvol_vars,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/HVOL/FINAL_HVOL_VARIANTS_FILTERED_COLLAPSED.tsv",row.names=F,col.names=T,sep="\t",quote=F)


####HCM####

#import HCM vars and remove duplicate entries
hcm_vars = read.table("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/HCM/FINAL_HCM_VARIANTS.tsv",header=T,sep="\t")
hcm_vars = unique(hcm_vars)

#keep only those above threshold
hcm_vars = subset(hcm_vars,hcm_vars$QUAL>=minqual)

#fix order of columns etc
hcm_vars$CHROM = gsub("chr","",hcm_vars$CHROM)
hcm_vars = cbind(hcm_vars[,c(1:6,8:9,7,10:12)],rep("",nrow(hcm_vars)),rep("",nrow(hcm_vars)),rep("",nrow(hcm_vars)),hcm_vars[,c(13:21)],rep("HCM",nrow(hcm_vars)),rep("GATKgcnv",nrow(hcm_vars)))
colnames(hcm_vars)[c(9,13:15,25,26)] = c("GATKgcnv_QUAL","AC","AC_Het","N_SAMPLES","COHORT","CALLER")

#using dplyr, collapse rows referring to the same variant in different samples
collapsed_hcm_vars = hcm_vars %>%
  group_by(across(-c(GATKgcnv_QUAL,SAMPLES))) %>%
  summarise(
    across(everything(), ~paste(., collapse = ",")),
    .groups = "drop"
  )
collapsed_hcm_vars = collapsed_hcm_vars[,colnames(hcm_vars)]

write.table(hcm_vars,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/HCM/FINAL_HCM_VARIANTS_FILTERED.tsv",row.names=F,col.names=T,sep="\t",quote=F)
write.table(collapsed_hcm_vars,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/HCM/FINAL_HCM_VARIANTS_FILTERED_COLLAPSED.tsv",row.names=F,col.names=T,sep="\t",quote=F)


####DCM####

#import DCM vars and remove duplicate entries
dcm_vars = read.table("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/DCM/FINAL_DCM_VARIANTS.tsv",header=T,sep="\t")
dcm_vars = unique(dcm_vars)

#keep only those above threshold
dcm_vars = subset(dcm_vars,dcm_vars$QUAL>=minqual)

#fix order of columns etc
dcm_vars$CHROM = gsub("chr","",dcm_vars$CHROM)
dcm_vars = cbind(dcm_vars[,c(1:6,8:9,7,10:12)],rep("",nrow(dcm_vars)),rep("",nrow(dcm_vars)),rep("",nrow(dcm_vars)),dcm_vars[,c(13:21)],rep("DCM",nrow(dcm_vars)),rep("GATKgcnv",nrow(dcm_vars)))
colnames(dcm_vars)[c(9,13:15,25,26)] = c("GATKgcnv_QUAL","AC","AC_Het","N_SAMPLES","COHORT","CALLER")

#using dplyr, collapse rows referring to the same variant in different samples
collapsed_dcm_vars = dcm_vars %>%
  group_by(across(-c(GATKgcnv_QUAL,SAMPLES))) %>%
  summarise(
    across(everything(), ~paste(., collapse = ",")),
    .groups = "drop"
  )
collapsed_dcm_vars = collapsed_dcm_vars[,colnames(dcm_vars)]

dcm_vars$CHROM = gsub("chr","",dcm_vars$CHROM)
dcm_vars = cbind(dcm_vars[,c(1:8,20,9:11)],rep("",nrow(dcm_vars)),rep("",nrow(dcm_vars)),rep("",nrow(dcm_vars)),dcm_vars[,c(21,12:19)])

write.table(dcm_vars,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/DCM/FINAL_DCM_VARIANTS_FILTERED.tsv",row.names=F,col.names=T,sep="\t",quote=F)
write.table(collapsed_dcm_vars,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/DCM/FINAL_DCM_VARIANTS_FILTERED_COLLAPSED.tsv",row.names=F,col.names=T,sep="\t",quote=F)



