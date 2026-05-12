library("ExomeDepth")

#define cohort to analyse
#########################
#note that the cohort name should be the name of the folder where case data are stored
cohort="DCM" #set to DCM, HCM or HVOL
#########################

#define genes of interest
geneset = c("ACTC1","ACTN2","BAG3","CSRP3","DES","DSP","JPH2","LMNA","MYBPC3","MYH7","MYL2","MYL3","NEXN","PLN","RBM20","SCN5A","TNNC1","TNNI3","TNNT2","TPM1","TRIM63","TTN","VCL")

#import target data and add exon identifiers
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

#define reference genome
fasta="/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta"

#generate counts objects and rename samples only in absence of a saved workspace
if(!paste("ExomeDepth_",cohort,".RData",sep="") %in% list.files("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts")){
  
  setwd(paste("/rds/general/user/fmazzaro/ephemeral/",cohort,sep=""))
  
  if(cohort == "DCM" | cohort == "HCM"){
    bamfiles_cases=subset(list.files(),!grepl("bai",list.files()))
    bamfiles_controls = paste("/rds/general/user/fmazzaro/ephemeral/HVOL/",subset(list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"),!grepl("bai",list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"))),sep="")
    bamfiles = c(bamfiles_cases,bamfiles_controls)
  }
  else{
    hvolfiles = subset(list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"),!grepl("bai",list.files("/rds/general/user/fmazzaro/ephemeral/HVOL")))
    toinclude = which(hvolfiles %in% c("14DB03367.bam","20SH04518.bam","14AD02262.bam","20QP07591.bam")) #these are hvol to be used as cases because other callers detected variants in them
    hvols_indexes = 1:1805
    hvols_left = subset(hvols_indexes,!hvols_indexes %in% toinclude)
    set.seed(13091998) #set seed for reproducibility
    hvols_picked_cases = sort(c(toinclude,sample(hvols_left,898,replace=F)))
    hvols_controls = subset(hvols_indexes,!hvols_indexes %in% hvols_picked_cases)
    bamfiles_cases=paste("/rds/general/user/fmazzaro/ephemeral/HVOL/",subset(list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"),!grepl("bai",list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"))),sep="")[hvols_picked_cases]
    bamfiles_controls = paste("/rds/general/user/fmazzaro/ephemeral/HVOL/",subset(list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"),!grepl("bai",list.files("/rds/general/user/fmazzaro/ephemeral/HVOL"))),sep="")[hvols_controls]
    bamfiles = c(bamfiles_cases,bamfiles_controls)
  }
  my.counts = getBamCounts(bed.frame=target_bed,bam.files=bamfiles,include.chr=FALSE,referenceFasta=fasta)
  new.colnames = substr(colnames(my.counts)[6:ncol(my.counts)],2,10)
  colnames(my.counts)[6:ncol(my.counts)] = new.colnames
  
  #clean the workspace by keeping only the my.counts object
  rm(list=setdiff(ls(),c("my.counts","bamfiles","bamfiles_cases","bamfiles_controls")))
  
  save.image(paste("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_",cohort,".RData",sep=""))
}else{
  load(paste("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_",cohort,".RData",sep=""))}

#remove chr prefix from the bed table
target_bed$chromosome = as.character(substr(as.character(target_bed$chromosome),4,5))

#select the most appropriate reference set from HVOL for each case and call CNVs
cnvlist = c()

for(i in 6:(6+length(bamfiles_cases)-1)){

  print(colnames(my.counts)[i])

  test.sample.name = colnames(my.counts)[i]
  ref.sample.names = colnames(my.counts)[(6+length(bamfiles_cases)):ncol(my.counts)]
  test.sample = as.vector(my.counts[,test.sample.name])
  ref.samples = as.matrix(my.counts[,ref.sample.names])

  #identify ideal reference samples
  my.choice = select.reference.set(test.counts=test.sample,reference.counts=ref.samples)

  if(length(my.choice$reference.choice)==1){
    my.matrix = as.matrix(my.counts[,my.choice$reference.choice,drop=FALSE])
  }
  else{
    my.matrix = as.matrix(my.counts[,my.choice$reference.choice])
  }
  my.reference.selected = apply(X=my.matrix,MAR=1,FUN=sum)

  #apply beta-binomial distribution to the full set of exons
  all.exons = new("ExomeDepth",test=test.sample,reference=my.reference.selected,formula="cbind(test, reference) ~ 1")

  #call CNVs and sort them by confidence level
  #all.exons = CallCNVs(x=all.exons,transition.probability=10^-4,chromosome=substr(as.character(my.counts$chromosome),4,5),start=my.counts$start,end=my.counts$end,name=my.counts$exon)
  all.exons = CallCNVs(x=all.exons,transition.probability=10^-4,chromosome=substr(as.character(my.counts$chromosome),4,5),start=my.counts$start,end=my.counts$end,name=my.counts$exon)
  
  if(length(all.exons@CNV.calls$BF)>0){
    cnv.calls = all.exons@CNV.calls[order(all.exons@CNV.calls$BF,decreasing=TRUE),]

    #annotate CNV calls manually (as the AnnotateExtra function from ExomeDepth seems to add NAs instead of gene names)
    sampleID = c()
    genes = c()
    exons = c()

    for(cnv in 1:nrow(cnv.calls)){
      sampleID = c(sampleID,test.sample.name)
      exonsinvolved = subset(target_bed,(target_bed$chromosome==cnv.calls$chromosome[cnv] & ((target_bed$start<=cnv.calls$start[cnv] & target_bed$end>=cnv.calls$start[cnv]) | (target_bed$start>=cnv.calls$start[cnv] & target_bed$end<=cnv.calls$end[cnv]) | (target_bed$start<=cnv.calls$start[cnv] & target_bed$end>=cnv.calls$end[cnv]) | (target_bed$start<=cnv.calls$start[cnv] & target_bed$end>=cnv.calls$end[cnv]))))
      exonids = paste(exonsinvolved$exon_id,collapse=",")
      genesinvolved = paste(unique(exonsinvolved$name),collapse=",")
      genes = c(genes,genesinvolved)
      exons = c(exons,exonids)

      #once processed the last CNV, keep only CNVs that hit genes of interest
      if(cnv == nrow(cnv.calls)){
        tokeep = c()
        for(y in 1:nrow(cnv.calls)){
          genesofcnv = unlist(strsplit(as.character(genes[y]),","))
          if(sum(genesofcnv %in% geneset)>0){
            tokeep = c(tokeep,TRUE)
          }
          else{
            tokeep = c(tokeep,FALSE)
          }
        }
      }
    }
    
    cnv.calls = cnv.calls[as.logical(tokeep),]
    sampleID = sampleID[as.logical(tokeep)]
    genes = genes[as.logical(tokeep)]
    exons = exons[as.logical(tokeep)]
    
    if(length(genesinvolved)>0){
      cnv.calls = cbind(cnv.calls,sampleID,genes,exons)
      cnvlist = rbind(cnvlist,cnv.calls)
    }
  }
}

#create a file with only high-PSI variants in TTN
write.table(cnvlist,paste("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_",cohort,"_l_RESULTS_ttnall.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
ttn_ex = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/ttnexons/TTN_highPSI_GRCh38.bed",sep="\t")
colnames(ttn_ex) = c("CHROM","START","END","EX_NO","PSI_DCM")
cnvlist$start = as.numeric(as.character(cnvlist$start))
cnvlist$end = as.numeric(as.character(cnvlist$end))
ttn_ex$START = as.numeric(as.character(ttn_ex$START))
ttn_ex$END = as.numeric(as.character(ttn_ex$END))
ttn_ex$PSI_DCM = as.numeric(as.character(ttn_ex$PSI_DCM))
  
edttn = subset(cnvlist,grepl("TTN",cnvlist$genes))
ednottn = subset(cnvlist,!grepl("TTN",cnvlist$genes))
  
tokeep = c()
for(v in 1:nrow(edttn)){
  n_ex = length(which((ttn_ex$START<=edttn$start[v] & ttn_ex$END>=edttn$end[v]) | (ttn_ex$START>=edttn$start[v] & ttn_ex$START<=edttn$end[v] & ttn_ex$END>=edttn$end[v]) |  (ttn_ex$END<=edttn$end[v] & ttn_ex$END>=edttn$start[v] & ttn_ex$START<=edttn$start[v]) | (ttn_ex$START>=edttn$start[v] & ttn_ex$END<=edttn$end[v])))
  if(n_ex>0){
    tokeep = c(tokeep,T)
  }
  else{
    tokeep = c(tokeep,F)
  }
}
edttn = edttn[as.logical(tokeep),]
cnvlist_ttnhigh = rbind(edttn,ednottn)
cnvlist_ttnhigh = cnvlist_ttnhigh[order(cnvlist_ttnhigh$chromosome,cnvlist_ttnhigh$start,decreasing=F),]
# cnvlist_ttnhigh = subset(cnvlist_ttnhigh,as.numeric(as.character(cnvlist_ttnhigh$BF))>4.472756)

write.table(cnvlist_ttnhigh,paste("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_",cohort,"_l_RESULTS_ttnhighpsi_expsize500.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


rm(list=setdiff(ls(),c("my.counts","bamfiles","bamfiles_cases","bamfiles_controls","cnvlist","cnvlist_ttnhigh","cohort")))
save.image(paste("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_",cohort,".RData",sep=""))
