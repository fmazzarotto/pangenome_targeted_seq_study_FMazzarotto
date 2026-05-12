library("ExomeDepth")

#define reference genome
fasta="/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta"

#define target regions
target_bed = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/overlap_AgiSSv5_HG002hc_AgiSSv8Sarcoid.bed")
gene_name = sapply(1:nrow(target_bed),function(x) unlist(strsplit(as.character(target_bed$V4[x]),"[,|]+"))[2])
target_bed = cbind(target_bed[,1:3],gene_name)

if(!"ExomeDepth_HG002_sarc.RData" %in% list.files("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts")){
  
  bam_hg002 = "/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/GIAB_SRA_WES_samples/HG002.bam"
  
  controlfiles = subset(list.files("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/wes_ref_controls/sarcoidosis/"),!grepl("bai",list.files("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/wes_ref_controls/sarcoidosis/")))
  controlfiles = paste0("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/data/wes_ref_controls/sarcoidosis/",controlfiles)
  
  bamfiles = c(bam_hg002,controlfiles)
  
  my.counts = getBamCounts(bed.frame=target_bed,bam.files=bamfiles,include.chr=F,referenceFasta=fasta)
  
  #select the most appropriate reference set from HVOL for HG002 and call CNVs
  test_sample = as.vector(my.counts[,6])
  controls = as.matrix(my.counts[,7:ncol(my.counts)])
  
  #clean the workspace by keeping only the my.counts object
  rm(list=setdiff(ls(),c("test_sample","my.choice","my.counts","bamfiles","bamfiles_cases","bamfiles_controls")))
  
  save.image("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_HG002_sarc.RData")
} else{
  load("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_HG002_sarc.RData")
}

#select the most appropriate reference control sample for HG002 and call CNVs
cnvlist = c()
test.sample.name = "HG002.bam"
if("exon" %in% colnames(my.counts)){
  ref.sample.names = colnames(my.counts)[7:ncol(my.counts)]
}else{
  ref.sample.names = colnames(my.counts)[6:ncol(my.counts)]
}

test.sample = as.vector(my.counts[,test.sample.name])
ref.samples = as.matrix(my.counts[,ref.sample.names])
my.choice = select.reference.set(test.counts=test.sample,reference.counts=ref.samples)


if(length(my.choice$reference.choice)==1){
  my.matrix = as.matrix(my.counts[,my.choice$reference.choice,drop=FALSE])
} else{
  my.matrix = as.matrix(my.counts[,my.choice$reference.choice])
}
my.reference.selected = apply(X=my.matrix,MAR=1,FUN=sum)

#apply beta-binomial distribution to the full set of exons
all.exons = new("ExomeDepth",test=test_sample,reference=my.reference.selected,formula="cbind(test, reference) ~ 1")

#call CNVs
all.exons = CallCNVs(x=all.exons,transition.probability=10^-4,chromosome=substr(as.character(my.counts$chromosome),4,5),start=my.counts$start,end=my.counts$end,name=rep(NA,nrow(my.counts)),expected.CNV.length = 50)
#all.exons = CallCNVs(x=all.exons,transition.probability=10^-4,chromosome=substr(as.character(my.counts$chromosome),4,5),start=my.counts$start,end=my.counts$end,name=rep(NA,nrow(my.counts)),expected.CNV.length = 500)
#all.exons = CallCNVs(x=all.exons,transition.probability=10^-4,chromosome=substr(as.character(my.counts$chromosome),4,5),start=my.counts$start,end=my.counts$end,name=rep(NA,nrow(my.counts)))

#save results
dir.create("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/exomedepth",showWarnings=F)
#write.table(all.exons@CNV.calls,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/HG002_ED_calls.tsv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(all.exons@CNV.calls,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/exomedepth/HG002_ED_calls_expsize50bp_sarcWES.tsv",col.names=T,row.names=F,sep="\t",quote=F)

#create and save BED file for intersection with truth set
bed = all.exons@CNV.calls[,c(7,5,6,3,9,12)]
bed$chromosome = paste0("chr",bed$chromosome)
bed$genotype=ifelse(bed$reads.ratio<0.25,"1/1",ifelse(bed$reads.ratio>=0.25 & bed$reads.ratio<0.75,"0/1",ifelse(bed$reads.ratio>=1.25 & bed$reads.ratio<1.75,"0/1",ifelse(bed$reads.ratio>=1.75,"1/1","0/0"))))
#write.table(bed,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/HG002_ED_calls.bed",col.names=F,row.names=F,sep="\t",quote=F)
write.table(bed,"/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/benchmarking_on_HG002/exomedepth/HG002_ED_calls_expsize50bp_sarcWES.bed",col.names=F,row.names=F,sep="\t",quote=F)
