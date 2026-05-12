
library(stringr)

genemap_file = "/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/genelists/genes.tsv"
genemap = read.table(genemap_file,header=T)

#set working directory
setwd("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/lumpy_results")

#modify VCF headers
system("zless HVOL_lumpy_raw.vcf.gz | sed s/#CHROM/CHROM/g > HVOL_temp.vcf")
system("zless DCM_lumpy_raw.vcf.gz | sed s/#CHROM/CHROM/g > DCM_temp.vcf")
system("zless HCM_lumpy_raw.vcf.gz | sed s/#CHROM/CHROM/g > HCM_temp.vcf")

#import tables
dcm = read.table("DCM_temp.vcf",header=T)
hcm = read.table("HCM_temp.vcf",header=T)
hvol = read.table("HVOL_temp.vcf",header=T)

#HVOL
#define new columns - CHROM
CHROM=hvol$CHROM
#START
START_hg38 = sapply(1:nrow(hvol),function(x) hvol$POS[x]+1)
#END
END_hg38 = as.numeric(as.character(gsub("END=","",sapply(1:nrow(hvol),function(x) strsplit(hvol$INFO[x],";")[[1]][3]))))
#REF and ALT
REF = rep(NA,nrow(hvol))
ALT = rep(NA,nrow(hvol))
#TYPE
TYPE = ifelse(hvol$ALT=="<DUP>","gain","loss")
#SVLEN
SVLEN = as.numeric(as.character(gsub("SVLEN=","",sapply(1:nrow(hvol),function(x) strsplit(hvol$INFO[x],";")[[1]][2]))))
#ED BF
ED_BF = rep(NA,nrow(hvol))
#POS_ORIG
POS_ORIG = hvol$POS
#REF_ORIG and ALT_ORIG
REF_ORIG = rep(NA,nrow(hvol))
ALT_ORIG = rep(NA,nrow(hvol))

#AC
count_nonref = function(gt) {
  # extract the first field before ":" (genotype string)
  g = sub(":.*", "", gt)
  # split into alleles
  alleles = unlist(strsplit(g, "[/|]"))
  # count how many are not "0"
  sum(alleles != "0")
}
#apply row-wise across genotype columns
AC = apply(hvol[,10:ncol(hvol)],1,function(row) {
  sum(sapply(row,count_nonref))
})

#AC_Het
AC_Het = rep(NA,nrow(hvol))

#N_SAMPLES
is_carrier = function(gt) {
  g = sub(":.*", "", gt)              # take genotype part before ":"
  alleles = unlist(strsplit(g, "[/|]"))
  # return 1 if any allele is non-reference (not 0, not missing)
  as.integer(any(alleles != "0"))
}

# Apply row-wise to columns 10–25
N_SAMPLES = apply(hvol[,10:ncol(hvol)],1,function(row) {
  sum(sapply(row, is_carrier))
})

#SAMPLES
get_carriers = function(row, sample_ids) {
  carriers = sapply(seq_along(row), function(i) {
    gt = sub(":.*", "", row[i])           # take genotype (before first ":")
    alleles = unlist(strsplit(gt, "[/|]"))
    if (any(alleles != "0" & alleles != ".")) {
      return(sample_ids[i])  # carrier
    } else {
      return(NA)             # not a carrier
    }
  })
  # Drop NAs and collapse into a single string (comma-separated IDs)
  paste(na.omit(carriers), collapse = ",")
}

# Apply to your dataframe
sample_ids = substr(colnames(hvol)[10:ncol(hvol)],2,10)

SAMPLES = apply(hvol[,10:ncol(hvol)], 1, function(row) {
  get_carriers(row, sample_ids)
})

#GENE,ENSG,ENST
GENE = c()
ENSG_ID = c()
ENST_ID = c()
for(x in 1:nrow(hvol)){
  gene = genemap[which(genemap$CHR == hvol$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"GENE"]
  GENE = c(GENE,gene)
  ensg = genemap[which(genemap$CHR == hvol$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"ENSG"]
  ENSG_ID = c(ENSG_ID,ensg)
  enst = genemap[which(genemap$CHR == hvol$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"ENST"]
  ENST_ID = c(ENST_ID,enst)
}

#CSQ, HVGSc, HGVSp
CSQ = rep(NA,nrow(hvol))
HGVSc = rep(NA,nrow(hvol))
HGVSp = rep(NA,nrow(hvol))

final_table = as.data.frame(cbind(CHROM,START_hg38,END_hg38,REF,ALT,TYPE,SVLEN,ED_BF,POS_ORIG,REF_ORIG,ALT_ORIG,AC,AC_Het,N_SAMPLES,SAMPLES,GENE,ENSG_ID,ENST_ID,CSQ,HGVSc,HGVSp))

#EXONS
#the add_exon_info function adds exons info to a table of variants. dataset=input table, objname=name of the resulting table.
add_exon_info <- function(dataset,objname){
  target_bed = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.bed")
  colnames(target_bed) = c("chromosome","start","end","name")
  exon_id = c()
  for(g in unique(target_bed$name)){
    ei = c()
    gs = subset(target_bed,target_bed$name==g)
    for(i in 1:nrow(gs)){
      ei = c(ei,paste(g,"_",i,sep=""))
    }
    exon_id = c(exon_id,ei)
  }
  target_bed = cbind(target_bed,exon_id)
  EXONS = c()
  for(cnv in 1:nrow(dataset)){
    exonsinvolved = subset(target_bed,(target_bed$chromosome==dataset$CHROM[cnv] & ((target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$START_hg38[cnv]) | (target_bed$start>=dataset$START_hg38[cnv] & target_bed$end<=dataset$END_hg38[cnv]) | (target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$END_hg38[cnv]) | (target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$END_hg38[cnv]))))
    exonids = paste(exonsinvolved$exon_id,collapse=",")
    EXONS = c(EXONS,exonids)
  }
  dataset = cbind(dataset,EXONS)
  
  #NOTE: this line removes variants not altering any exon
  dataset = subset(dataset,dataset$EXONS != "")
  
  assign(objname,dataset,envir=.GlobalEnv)
}

add_exon_info(final_table,"final_table")

#transform EXON list in a range
exon_range = c()
temp_list = sapply(1:nrow(final_table),function(x) strsplit(final_table$EXONS[x],","))

for(l in 1:length(temp_list)){
  min_exon = min(as.numeric(gsub(".*_","",temp_list[[l]])))
  max_exon = max(as.numeric(gsub(".*_","",temp_list[[l]])))
  
  if(min_exon==max_exon){
    exon_range = c(exon_range,min_exon)
  }
  else{
    exon_range = c(exon_range,paste0(min_exon,"-",max_exon))
  }
}
final_table$EXONS = exon_range

#gnomAD_AF
gnomad = read.csv("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gnomad/gnomAD_SV_CMgenes.csv",header=T,colClasses = c(rep("numeric",5),rep("character",3),rep("numeric",6)))
gnomad$CHROM=paste0("chr",gnomad$CHROM)
if(!"EXONS" %in% colnames(gnomad)){
  add_exon_info(gnomad,"gnomad")
}
gnomad_common = subset(gnomad,gnomad$Allele.Frequency >= 0.0001)

gnomAD_AF = c()
for(v in 1:nrow(final_table)){
  
  if(sum(gnomad$EXONS == final_table$EXONS[v]) > 0){
    gnomAD_AF = c(gnomAD_AF,"Check_manually")
  }
  else{
    gnomAD_AF = c(gnomAD_AF,0)
  }
}
final_table = cbind(final_table,gnomAD_AF)

#COHORT
COHORT = rep("HVOL",nrow(final_table))
final_table = cbind(final_table,COHORT)

#CALLER
COHORT = rep("LUMPY-SV",nrow(final_table))

#CENTRE
CENTRE=c()
sample_list = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/sample_lists/London_Aswan_sample_list_complete.tsv",header=T)
for(x in 1:nrow(final_table)){
  if(final_table$N_SAMPLES[x] > 1){
    CENTRE=c(CENTRE,"Check manually")
  }
  else{
    if(startsWith(final_table$SAMPLES[x],"CV")){
      CENTRE=c(CENTRE,"Florence")
    }
    else if(grepl("FNM",final_table$SAMPLES[x]) | startsWith(final_table$SAMPLES[x],"KM")){
      CENTRE=c(CENTRE,"Prague")
    }
    else{
      CENTRE=c(CENTRE,sample_list[which(sample_list$BruNumber==final_table$SAMPLES[x]),"CENTRE"])
    }
  }
}

final_table = cbind(final_table,CENTRE)
final_table$CHROM = gsub("chr","",final_table$CHROM)
hvol_table = final_table
rm(final_table)

###############################################################################################################

#DCM
#define new columns - CHROM
CHROM=dcm$CHROM
#START
START_hg38 = sapply(1:nrow(dcm),function(x) dcm$POS[x]+1)
#END
END_hg38 = as.numeric(as.character(gsub("END=","",sapply(1:nrow(dcm),function(x) strsplit(dcm$INFO[x],";")[[1]][3]))))
#REF and ALT
REF = rep(NA,nrow(dcm))
ALT = rep(NA,nrow(dcm))
#TYPE
TYPE = ifelse(dcm$ALT=="<DUP>","gain","loss")
#SVLEN
SVLEN = as.numeric(as.character(gsub("SVLEN=","",sapply(1:nrow(dcm),function(x) strsplit(dcm$INFO[x],";")[[1]][2]))))
#ED BF
ED_BF = rep(NA,nrow(dcm))
#POS_ORIG
POS_ORIG = dcm$POS
#REF_ORIG and ALT_ORIG
REF_ORIG = rep(NA,nrow(dcm))
ALT_ORIG = rep(NA,nrow(dcm))

#AC
count_nonref = function(gt) {
  # extract the first field before ":" (genotype string)
  g = sub(":.*", "", gt)
  # split into alleles
  alleles = unlist(strsplit(g, "[/|]"))
  # count how many are not "0"
  sum(alleles != "0")
}
#apply row-wise across genotype columns
AC = apply(dcm[,10:ncol(dcm)],1,function(row) {
  sum(sapply(row,count_nonref))
})

#AC_Het
AC_Het = rep(NA,nrow(dcm))

#N_SAMPLES
is_carrier = function(gt) {
  g = sub(":.*", "", gt)              # take genotype part before ":"
  alleles = unlist(strsplit(g, "[/|]"))
  # return 1 if any allele is non-reference (not 0, not missing)
  as.integer(any(alleles != "0"))
}

# Apply row-wise to columns 10–25
N_SAMPLES = apply(dcm[,10:ncol(dcm)],1,function(row) {
  sum(sapply(row, is_carrier))
})

#SAMPLES
get_carriers = function(row, sample_ids) {
  carriers = sapply(seq_along(row), function(i) {
    gt = sub(":.*", "", row[i])           # take genotype (before first ":")
    alleles = unlist(strsplit(gt, "[/|]"))
    if (any(alleles != "0" & alleles != ".")) {
      return(sample_ids[i])  # carrier
    } else {
      return(NA)             # not a carrier
    }
  })
  # Drop NAs and collapse into a single string (comma-separated IDs)
  paste(na.omit(carriers), collapse = ",")
}

# Apply to your dataframe
sample_ids = substr(colnames(dcm)[10:ncol(dcm)],2,10)

SAMPLES = apply(dcm[,10:ncol(dcm)], 1, function(row) {
  get_carriers(row, sample_ids)
})

#GENE,ENSG,ENST
GENE = c()
ENSG_ID = c()
ENST_ID = c()
for(x in 1:nrow(dcm)){
  gene = genemap[which(genemap$CHR == dcm$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"GENE"]
  GENE = c(GENE,gene)
  ensg = genemap[which(genemap$CHR == dcm$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"ENSG"]
  ENSG_ID = c(ENSG_ID,ensg)
  enst = genemap[which(genemap$CHR == dcm$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"ENST"]
  ENST_ID = c(ENST_ID,enst)
}

#CSQ, HVGSc, HGVSp
CSQ = rep(NA,nrow(dcm))
HGVSc = rep(NA,nrow(dcm))
HGVSp = rep(NA,nrow(dcm))

final_table = as.data.frame(cbind(CHROM,START_hg38,END_hg38,REF,ALT,TYPE,SVLEN,ED_BF,POS_ORIG,REF_ORIG,ALT_ORIG,AC,AC_Het,N_SAMPLES,SAMPLES,GENE,ENSG_ID,ENST_ID,CSQ,HGVSc,HGVSp))

#EXONS
#the add_exon_info function adds exons info to a table of variants. dataset=input table, objname=name of the resulting table.
add_exon_info <- function(dataset,objname){
  target_bed = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.bed")
  colnames(target_bed) = c("chromosome","start","end","name")
  exon_id = c()
  for(g in unique(target_bed$name)){
    ei = c()
    gs = subset(target_bed,target_bed$name==g)
    for(i in 1:nrow(gs)){
      ei = c(ei,paste(g,"_",i,sep=""))
    }
    exon_id = c(exon_id,ei)
  }
  target_bed = cbind(target_bed,exon_id)
  EXONS = c()
  for(cnv in 1:nrow(dataset)){
    exonsinvolved = subset(target_bed,(target_bed$chromosome==dataset$CHROM[cnv] & ((target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$START_hg38[cnv]) | (target_bed$start>=dataset$START_hg38[cnv] & target_bed$end<=dataset$END_hg38[cnv]) | (target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$END_hg38[cnv]) | (target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$END_hg38[cnv]))))
    exonids = paste(exonsinvolved$exon_id,collapse=",")
    EXONS = c(EXONS,exonids)
  }
  dataset = cbind(dataset,EXONS)
  
  #NOTE: this line removes variants not altering any exon
  dataset = subset(dataset,dataset$EXONS != "")
  
  assign(objname,dataset,envir=.GlobalEnv)
}

add_exon_info(final_table,"final_table")

#transform EXON list in a range
exon_range = c()
temp_list = sapply(1:nrow(final_table),function(x) strsplit(final_table$EXONS[x],","))

for(l in 1:length(temp_list)){
  min_exon = min(as.numeric(gsub(".*_","",temp_list[[l]])))
  max_exon = max(as.numeric(gsub(".*_","",temp_list[[l]])))
  
  if(min_exon==max_exon){
    exon_range = c(exon_range,min_exon)
  }
  else{
    exon_range = c(exon_range,paste0(min_exon,"-",max_exon))
  }
}
final_table$EXONS = exon_range

#gnomAD_AF
gnomad = read.csv("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gnomad/gnomAD_SV_CMgenes.csv",header=T,colClasses = c(rep("numeric",5),rep("character",3),rep("numeric",6)))
gnomad$CHROM=paste0("chr",gnomad$CHROM)
if(!"EXONS" %in% colnames(gnomad)){
  add_exon_info(gnomad,"gnomad")
}
gnomad_common = subset(gnomad,gnomad$Allele.Frequency >= 0.0001)

gnomAD_AF = c()
for(v in 1:nrow(final_table)){
  
  if(sum(gnomad$EXONS == final_table$EXONS[v]) > 0){
    gnomAD_AF = c(gnomAD_AF,"Check_manually")
  }
  else{
    gnomAD_AF = c(gnomAD_AF,0)
  }
}
final_table = cbind(final_table,gnomAD_AF)

#COHORT
COHORT = rep("DCM",nrow(final_table))
final_table = cbind(final_table,COHORT)

#CALLER
COHORT = rep("LUMPY-SV",nrow(final_table))

#CENTRE
CENTRE=c()
sample_list = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/sample_lists/London_Aswan_sample_list_complete.tsv",header=T)
for(x in 1:nrow(final_table)){
  if(final_table$N_SAMPLES[x] > 1){
    CENTRE=c(CENTRE,"Check manually")
  }
  else{
    if(startsWith(final_table$SAMPLES[x],"CV")){
      CENTRE=c(CENTRE,"Florence")
    }
    else if(grepl("FNM",final_table$SAMPLES[x]) | startsWith(final_table$SAMPLES[x],"KM")){
      CENTRE=c(CENTRE,"Prague")
    }
    else{
      CENTRE=c(CENTRE,sample_list[which(sample_list$BruNumber==final_table$SAMPLES[x]),"CENTRE"])
    }
  }
}

final_table = cbind(final_table,CENTRE)
final_table$CHROM = gsub("chr","",final_table$CHROM)
dcm_table = final_table
rm(final_table)

###############################################################################################################

#HCM
#define new columns - CHROM
CHROM=hcm$CHROM
#START
START_hg38 = sapply(1:nrow(hcm),function(x) hcm$POS[x]+1)
#END
END_hg38 = as.numeric(as.character(gsub("END=","",sapply(1:nrow(hcm),function(x) strsplit(hcm$INFO[x],";")[[1]][3]))))
#REF and ALT
REF = rep(NA,nrow(hcm))
ALT = rep(NA,nrow(hcm))
#TYPE
TYPE = ifelse(hcm$ALT=="<DUP>","gain","loss")
#SVLEN
SVLEN = as.numeric(as.character(gsub("SVLEN=","",sapply(1:nrow(hcm),function(x) strsplit(hcm$INFO[x],";")[[1]][2]))))
#ED BF
ED_BF = rep(NA,nrow(hcm))
#POS_ORIG
POS_ORIG = hcm$POS
#REF_ORIG and ALT_ORIG
REF_ORIG = rep(NA,nrow(hcm))
ALT_ORIG = rep(NA,nrow(hcm))

#AC
count_nonref = function(gt) {
  # extract the first field before ":" (genotype string)
  g = sub(":.*", "", gt)
  # split into alleles
  alleles = unlist(strsplit(g, "[/|]"))
  # count how many are not "0"
  sum(alleles != "0")
}
#apply row-wise across genotype columns
AC = apply(hcm[,10:ncol(hcm)],1,function(row) {
  sum(sapply(row,count_nonref))
})

#AC_Het
AC_Het = rep(NA,nrow(hcm))

#N_SAMPLES
is_carrier = function(gt) {
  g = sub(":.*", "", gt)              # take genotype part before ":"
  alleles = unlist(strsplit(g, "[/|]"))
  # return 1 if any allele is non-reference (not 0, not missing)
  as.integer(any(alleles != "0"))
}

# Apply row-wise to columns 10–25
N_SAMPLES = apply(hcm[,10:ncol(hcm)],1,function(row) {
  sum(sapply(row, is_carrier))
})

#SAMPLES
get_carriers = function(row, sample_ids) {
  carriers = sapply(seq_along(row), function(i) {
    gt = sub(":.*", "", row[i])           # take genotype (before first ":")
    alleles = unlist(strsplit(gt, "[/|]"))
    if (any(alleles != "0" & alleles != ".")) {
      return(sample_ids[i])  # carrier
    } else {
      return(NA)             # not a carrier
    }
  })
  # Drop NAs and collapse into a single string (comma-separated IDs)
  paste(na.omit(carriers), collapse = ",")
}

# Apply to your dataframe
sample_ids = substr(colnames(hcm)[10:ncol(hcm)],2,10)

SAMPLES = apply(hcm[,10:ncol(hcm)], 1, function(row) {
  get_carriers(row, sample_ids)
})

#GENE,ENSG,ENST
GENE = c()
ENSG_ID = c()
ENST_ID = c()
for(x in 1:nrow(hcm)){
  gene = genemap[which(genemap$CHR == hcm$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"GENE"]
  GENE = c(GENE,gene)
  ensg = genemap[which(genemap$CHR == hcm$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"ENSG"]
  ENSG_ID = c(ENSG_ID,ensg)
  enst = genemap[which(genemap$CHR == hcm$CHROM[x] & ((genemap$START<START_hg38[x] & genemap$END>START_hg38[x]) | (genemap$START<END_hg38[x] & genemap$END>END_hg38[x]))),"ENST"]
  ENST_ID = c(ENST_ID,enst)
}

#CSQ, HVGSc, HGVSp
CSQ = rep(NA,nrow(hcm))
HGVSc = rep(NA,nrow(hcm))
HGVSp = rep(NA,nrow(hcm))

final_table = as.data.frame(cbind(CHROM,START_hg38,END_hg38,REF,ALT,TYPE,SVLEN,ED_BF,POS_ORIG,REF_ORIG,ALT_ORIG,AC,AC_Het,N_SAMPLES,SAMPLES,GENE,ENSG_ID,ENST_ID,CSQ,HGVSc,HGVSp))

#EXONS
#the add_exon_info function adds exons info to a table of variants. dataset=input table, objname=name of the resulting table.
add_exon_info <- function(dataset,objname){
  target_bed = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.bed")
  colnames(target_bed) = c("chromosome","start","end","name")
  exon_id = c()
  for(g in unique(target_bed$name)){
    ei = c()
    gs = subset(target_bed,target_bed$name==g)
    for(i in 1:nrow(gs)){
      ei = c(ei,paste(g,"_",i,sep=""))
    }
    exon_id = c(exon_id,ei)
  }
  target_bed = cbind(target_bed,exon_id)
  EXONS = c()
  for(cnv in 1:nrow(dataset)){
    exonsinvolved = subset(target_bed,(target_bed$chromosome==dataset$CHROM[cnv] & ((target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$START_hg38[cnv]) | (target_bed$start>=dataset$START_hg38[cnv] & target_bed$end<=dataset$END_hg38[cnv]) | (target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$END_hg38[cnv]) | (target_bed$start<=dataset$START_hg38[cnv] & target_bed$end>=dataset$END_hg38[cnv]))))
    exonids = paste(exonsinvolved$exon_id,collapse=",")
    EXONS = c(EXONS,exonids)
  }
  dataset = cbind(dataset,EXONS)
  
  #NOTE: this line removes variants not altering any exon
  dataset = subset(dataset,dataset$EXONS != "")
  
  assign(objname,dataset,envir=.GlobalEnv)
}

add_exon_info(final_table,"final_table")

#transform EXON list in a range
exon_range = c()
temp_list = sapply(1:nrow(final_table),function(x) strsplit(final_table$EXONS[x],","))

for(l in 1:length(temp_list)){
  min_exon = min(as.numeric(gsub(".*_","",temp_list[[l]])))
  max_exon = max(as.numeric(gsub(".*_","",temp_list[[l]])))
  
  if(min_exon==max_exon){
    exon_range = c(exon_range,min_exon)
  }
  else{
    exon_range = c(exon_range,paste0(min_exon,"-",max_exon))
  }
}
final_table$EXONS = exon_range

#gnomAD_AF
gnomad = read.csv("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gnomad/gnomAD_SV_CMgenes.csv",header=T,colClasses = c(rep("numeric",5),rep("character",3),rep("numeric",6)))
gnomad$CHROM=paste0("chr",gnomad$CHROM)
if(!"EXONS" %in% colnames(gnomad)){
  add_exon_info(gnomad,"gnomad")
}
gnomad_common = subset(gnomad,gnomad$Allele.Frequency >= 0.0001)

gnomAD_AF = c()
for(v in 1:nrow(final_table)){
  
  if(sum(gnomad$EXONS == final_table$EXONS[v]) > 0){
    gnomAD_AF = c(gnomAD_AF,"Check_manually")
  }
  else{
    gnomAD_AF = c(gnomAD_AF,0)
  }
}
final_table = cbind(final_table,gnomAD_AF)

#COHORT
COHORT = rep("HCM",nrow(final_table))
final_table = cbind(final_table,COHORT)

#CALLER
COHORT = rep("LUMPY-SV",nrow(final_table))

#CENTRE
CENTRE=c()
sample_list = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/sample_lists/London_Aswan_sample_list_complete.tsv",header=T)
for(x in 1:nrow(final_table)){
  if(final_table$N_SAMPLES[x] > 1){
    CENTRE=c(CENTRE,"Check manually")
  }
  else{
    if(startsWith(final_table$SAMPLES[x],"CV")){
      CENTRE=c(CENTRE,"Florence")
    }
    else if(grepl("FNM",final_table$SAMPLES[x]) | startsWith(final_table$SAMPLES[x],"KM")){
      CENTRE=c(CENTRE,"Prague")
    }
    else{
      CENTRE=c(CENTRE,sample_list[which(sample_list$BruNumber==final_table$SAMPLES[x]),"CENTRE"])
    }
  }
}

final_table = cbind(final_table,CENTRE)
final_table$CHROM = gsub("chr","",final_table$CHROM)
hcm_table = final_table

###############################################################################################################

final_table = rbind(hvol_table,dcm_table,hcm_table)
final_table = subset(final_table,abs(as.numeric(as.character(final_table$SVLEN)))>=20)
write.table(final_table,"lumpy_candidates_FINAL_pre_manual_edits.tsv",col.names=T,row.names=F,quote=F,sep="\t")

