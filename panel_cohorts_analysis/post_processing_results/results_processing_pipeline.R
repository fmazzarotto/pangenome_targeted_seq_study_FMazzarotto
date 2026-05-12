#import files (manta hcm and hvol files have got no variants and will not be imported - see hashed lines below)
gnomad = read.csv("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gnomad/gnomAD_SV_CMgenes.csv",header=T,colClasses = c(rep("numeric",5),rep("character",3),rep("numeric",6)))
ed_hcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_HCM_l_RESULTS_ttnhighpsi.tsv",header=T,colClasses = c(rep("numeric",2),"character",rep("numeric",4),"character",rep("numeric",4),rep("character",3)))
ed_dcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_DCM_l_RESULTS_ttnhighpsi.tsv",header=T,colClasses = c(rep("numeric",2),"character",rep("numeric",4),"character",rep("numeric",4),rep("character",3)))
ed_hvol = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_HVOL_l_RESULTS_ttnhighpsi.tsv",header=T,colClasses = c(rep("numeric",2),"character",rep("numeric",4),"character",rep("numeric",4),rep("character",3)))
gatk_dcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/gatk_dcm_l_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))
gatk_hcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/gatk_hcm_l_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))
gatk_hvol = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/gatk_hvol_l_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))
manta_dcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/manta_dcm_l_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",4),rep("character",7),"numeric"))
#manta_hcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/manta_hcm_l_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",4),rep("character",7),"numeric"))
#manta_hvol = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/manta_hvol_l_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",4),rep("character",7),"numeric"))
graf_dcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/graf_dcm_d_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))
graf_hcm = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/graf_hcm_h_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))
graf_hvol_d = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/graf_hvol_d_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))
graf_hvol_h = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/graf_hvol_h_RESULTS.tsv",header=T,colClasses = c("character","numeric",rep("character",3),rep("numeric",3),rep("character",7),"numeric"))

#define gene sets
dcmgeneset = c("BAG3","DES","LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN","DSP","ACTC1","ACTN2","JPH2","NEXN","TNNI3","TPM1","VCL")
hcmgeneset = c("MYBPC3","MYH7","TNNT2","TNNI3","MYL2","MYL3","ACTC1","TPM1","ACTN2","CSRP3","TNNC1","PLN","JPH2","TRIM63")


#####################################
######define relevant functions######
#####################################

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

#the add_start_end functions convert the POS notation to START and END, accordingly replacing the unchanged ref base to '-' (for compatibility with ExomeDepth START and END)
add_start_end_gatk_graf <- function(dataset,objname){
  START_hg38 = c()
  END_hg38 = c()
  TYPE = c()
  SVLEN = c()
  REF = c()
  ALT = c()
  colnames(dataset)[which(colnames(dataset)=="POS")] = "POS_orig"
  colnames(dataset)[which(colnames(dataset)=="REF")] = "REF_orig"
  colnames(dataset)[which(colnames(dataset)=="ALT")] = "ALT_orig"
  for(var in 1:nrow(dataset)){
    reflen = nchar(dataset$REF[var])
    altlen = nchar(dataset$ALT[var])
    if(reflen > altlen){
      if(altlen == 1){
        REF = c(REF,substr(dataset$REF[var],2,reflen))
        ALT = c(ALT,"-")
      }
      else{
        print(dataset[var,])
        stop("Check dataset: there are complex variants in the dataset for which the add_start_end function may not work correctly")
      }
      TYPE = c(TYPE,"loss")
      START_hg38 = c(START_hg38,dataset$POS[var]+1)
      END_hg38 = c(END_hg38,dataset$POS[var]+nchar(dataset$REF[var])-1)
      SVLEN = c(SVLEN,as.numeric(paste("-",nchar(REF[length(REF)]),sep="")))
    }
    else{
      if(reflen == 1){
        REF = c(REF,"-")
        ALT = c(ALT,substr(dataset$ALT[var],2,altlen))
      }
      else{
        stop("Check dataset: there are complex variants in the dataset for which the add_start_end function may not work correctly")
      }
      TYPE = c(TYPE,"gain")
      START_hg38 = c(START_hg38,dataset$POS[var])
      END_hg38 = c(END_hg38,dataset$POS[var]+1)
      SVLEN = c(SVLEN,nchar(ALT[length(ALT)]))
    }
  }
  
  BF = rep(NA,nrow(dataset))
  dataset = cbind(dataset,START_hg38,END_hg38,REF,ALT,TYPE,SVLEN,BF)
  dataset = dataset[,c(1,17:23,2:4,6:16)]
  colnames(dataset)[which(colnames(dataset) == "NS")] = "N_SAMPLES"
  colnames(dataset)[which(colnames(dataset) == "SYMBOL")] = "GENE"
  colnames(dataset)[which(colnames(dataset) == "Gene")] = "ENSG_ID"
  colnames(dataset)[which(colnames(dataset) == "Feature")] = "ENST_ID"
  colnames(dataset)[which(colnames(dataset) == "Consequence")] = "CSQ"

  assign(objname,dataset,envir=.GlobalEnv)
}

add_start_end_manta <- function(dataset,objname){
  START_hg38 = c()
  END_hg38 = c()
  TYPE = c()
  REF = c()
  ALT = c()
  colnames(dataset)[which(colnames(dataset)=="POS")] = "POS_orig"
  colnames(dataset)[which(colnames(dataset)=="REF")] = "REF_orig"
  colnames(dataset)[which(colnames(dataset)=="ALT")] = "ALT_orig"
  for(var in 1:nrow(dataset)){
    
    if(!grepl("<",dataset$ALT[var]) && !grepl("<",dataset$REF[var])){  #if ref and alt sequences are provided
      reflen = nchar(dataset$REF[var])
      altlen = nchar(dataset$ALT[var])
      
      if(substr(dataset$SVLEN[var],1,1) == "-"){
        if(altlen == 1){
          ALT = c(ALT,"-")
          REF = c(REF,substr(dataset$REF[var],2,reflen))
        }
        else{
          print(dataset[var,])
          stop("Check dataset: there are complex variants in the dataset for which the add_start_end function may not work correctly")
        }
        TYPE = c(TYPE,"loss")
        START_hg38 = c(START_hg38,dataset$POS[var]+1)
        END_hg38 = c(END_hg38,dataset$POS[var]+nchar(dataset$REF[var])-1)
      }
      else{
        if(reflen == 1){
          REF = c(REF,"-")
          ALT = c(ALT,substr(dataset$ALT[var],2,altlen))
        }
        else{
          stop("Check dataset: there are complex variants in the dataset for which the add_start_end function may not work correctly")
        }
        TYPE = c(TYPE,"gain")
        START_hg38 = c(START_hg38,dataset$POS[var])
        END_hg38 = c(END_hg38,dataset$POS[var]+1)
      }
    }
    
    else{
      if(substr(dataset$SVLEN[var],1,1) == "-"){ #if it is a deletion
        ALT = c(ALT,NA)
        REF = c(REF,NA)
        START_hg38 = c(START_hg38,dataset$POS[var]+1)
        END_hg38 = c(END_hg38,START_hg38[length(START_hg38)]+abs(dataset$SVLEN[var])-1)
        TYPE = c(TYPE,"loss")
      }
      else if(grepl("DUP",dataset$ALT[var])){ #if it is a duplication 
        ALT = c(ALT,NA)
        REF = c(REF,NA)
        START_hg38 = c(START_hg38,dataset$POS[var])
        END_hg38 = c(END_hg38,dataset$POS[var]+dataset$SVLEN[var])
        TYPE = c(TYPE,"gain")
      }
      else{ #if it is an insertion
        ALT = c(ALT,NA)
        REF = c(REF,NA)
        START_hg38 = c(START_hg38,dataset$POS[var])
        END_hg38 = c(END_hg38,dataset$POS[var]+1)
        TYPE = c(TYPE,"gain")
      }
    }
  }
  
  BF = rep(NA,nrow(dataset))
  dataset = cbind(dataset,START_hg38,END_hg38,REF,ALT,TYPE,BF)
  dataset = dataset[,c(1,18:22,6,23,2:4,7:17)]
  colnames(dataset)[which(colnames(dataset) == "NS")] = "N_SAMPLES"
  colnames(dataset)[which(colnames(dataset) == "SYMBOL")] = "GENE"
  colnames(dataset)[which(colnames(dataset) == "Gene")] = "ENSG_ID"
  colnames(dataset)[which(colnames(dataset) == "Feature")] = "ENST_ID"
  colnames(dataset)[which(colnames(dataset) == "Consequence")] = "CSQ"
  
  assign(objname,dataset,envir=.GlobalEnv)
}

#fix_ed_results_format converts the multiple rows of a recurrent ExomeDepth-called variant into a single row with all carriers listed. The median BF and read counts across carriers are displayed.
#it also adds missing info, renames columns according to colnames in other callers' datasets and keeps columns of interest 

fix_ed_results_format <- function(ed_dataset,objname){

  genes = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/genelists/genes.tsv",header=T)

  recur_vars = ed_dataset[as.logical((duplicated(ed_dataset$type) | duplicated(ed_dataset$type,fromLast = T)) & (duplicated(ed_dataset$id) | duplicated(ed_dataset$id,fromLast = T))),]
  nonrecur_vars = ed_dataset[!as.logical((duplicated(ed_dataset$type) | duplicated(ed_dataset$type,fromLast = T)) & (duplicated(ed_dataset$id) | duplicated(ed_dataset$id,fromLast = T))),]
  
  if(nrow(recur_vars) > 0){
    newlines = c()
    for(recurvar in unique(recur_vars$id)){
      varsubset = subset(recur_vars,recur_vars$id == recurvar)
      BF = median(varsubset$BF)
      sampleID = paste(varsubset$sampleID,collapse=",")
      reads.expected = median(varsubset$reads.expected)
      reads.observed = median(varsubset$reads.observed)
      reads.ratio = median(varsubset$reads.ratio)
      newline = c(varsubset$start.p[1],varsubset$end.p[1],varsubset$type[1],varsubset$nexons[1],varsubset$start[1],varsubset$end[1],varsubset$chromosome[1],varsubset$id[1],BF,reads.expected,reads.observed,reads.ratio,sampleID,varsubset$genes[1],varsubset$exons[1])
      newlines = rbind(newlines,newline)
    }
    newlines = as.data.frame(newlines)
    colnames(newlines) = colnames(nonrecur_vars)
    ed_dataset = rbind(nonrecur_vars,newlines)
    ed_dataset$chromosome = as.numeric(ed_dataset$chromosome)
    ed_dataset$start = as.numeric(ed_dataset$start)
    ed_dataset = ed_dataset[order(ed_dataset$chromosome,ed_dataset$start),]
    rownames(ed_dataset) = 1:nrow(ed_dataset)
  }
  
  REF = rep(NA,nrow(ed_dataset))
  ALT = rep(NA,nrow(ed_dataset))
  POS_orig = rep(NA,nrow(ed_dataset))
  REF_orig = rep(NA,nrow(ed_dataset))
  ALT_orig = rep(NA,nrow(ed_dataset))
  AC_Het = rep(NA,nrow(ed_dataset))
  AC = sapply(1:nrow(ed_dataset),function(v) length(unlist(strsplit(ed_dataset$sampleID[v],","))))
  NS = AC
  SVLEN = sapply(1:nrow(ed_dataset),function(v) as.numeric(ed_dataset$end[v])-as.numeric(ed_dataset$start[v])+1)
  Consequence = rep(NA,nrow(ed_dataset))
  HGVSp = rep(NA,nrow(ed_dataset))
  HGVSc = rep(NA,nrow(ed_dataset))
  
  if(!"EXONS" %in% colnames(gnomad)){
    gnomad = read.csv("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gnomad/gnomAD_SV_CMgenes.csv",header=T,colClasses = c(rep("numeric",5),rep("character",3),rep("numeric",6)))
    add_exon_info(gnomad,"gnomad")
  }
  gnomad_common = subset(gnomad,gnomad$Allele.Frequency >= 0.0001)
  
  Gene = c()
  Feature = c()
  gnomAD_AF = c()
  for(v in 1:nrow(ed_dataset)){
    
    if(sum(gnomad$EXONS == ed_dataset$exons[v]) > 0){
      gnomAD_AF = c(gnomAD_AF,"Check_manually")
    }
    else{
      gnomAD_AF = c(gnomAD_AF,0)
    }
    
    vargenes = unlist(strsplit(ed_dataset$genes,","))
    varensg = c()
    varenst = c()
    for(g in vargenes){
      varensg = c(varensg,genes[which(genes$GENE == g),"ENSG"])
      varenst = c(varensg,genes[which(genes$GENE == g),"ENST"])
    }
    varensg = paste(varensg,collapse=",")
    varenst = paste(varenst,collapse=",")
    Gene = c(Gene,varensg)
    Feature = c(Feature,varenst)
  }
  
  ed_dataset = cbind(ed_dataset,REF,ALT,POS_orig,REF_orig,ALT_orig,AC_Het,AC,NS,SVLEN,Consequence,HGVSc,HGVSp,Gene,Feature,gnomAD_AF)
  
  if("exons" %in% colnames(ed_dataset)){
    ed_dataset = ed_dataset[,c(7,5:6,16:17,3,24,9,18:20,22,21,23,13:15,28:29,25:27,30)]
  }
  else{
    ed_dataset = ed_dataset[,c(7,5:6,16:17,3,24,9,18:20,22,21,23,13:14,28:29,25:27,30)]
  }
  
  colnames(ed_dataset)[which(colnames(ed_dataset) == "chromosome")] = "CHROM"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "start")] = "START_hg38"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "end")] = "END_hg38"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "NS")] = "N_SAMPLES"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "sampleID")] = "SAMPLES"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "genes")] = "GENE"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "Gene")] = "ENSG_ID"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "Feature")] = "ENST_ID"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "exons")] = "EXONS"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "Consequence")] = "CSQ"
  colnames(ed_dataset)[which(colnames(ed_dataset) == "type")] = "TYPE"
  
  assign(objname,ed_dataset,envir=.GlobalEnv)
}

#add_column just adds a column filled with the same information for every variant in case needed, with a name chosen by the user
add_column <- function(dataset,objname,column_name,content){
  newcol = rep(content,nrow(dataset))
  dataset = cbind(dataset,newcol)
  colnames(dataset)[which(colnames(dataset) == "newcol")] = column_name
  
  assign(objname,dataset,envir=.GlobalEnv)
  }

#################################################################
######remove low-psi exon TTN variants from the gnomAD table#####
#################################################################
ttn_ex = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/ttnexons/TTN_highPSI_GRCh38.bed",sep="\t")
colnames(ttn_ex) = c("CHROM","START","END","EX_NO","PSI_DCM")
gnomad$START_hg38 = as.numeric(as.character(gnomad$START_hg38))
gnomad$END_hg38 = as.numeric(as.character(gnomad$END_hg38))
ttn_ex$START = as.numeric(as.character(ttn_ex$START))
ttn_ex$END = as.numeric(as.character(ttn_ex$END))
ttn_ex$PSI_DCM = as.numeric(as.character(ttn_ex$PSI_DCM))

gnomadttn = subset(gnomad,gnomad$GENE == "TTN")
gnomadnottn = subset(gnomad,gnomad$GENE != "TTN")

tokeep = c()
for(v in 1:nrow(gnomadttn)){
  n_ex = length(which((ttn_ex$START<=gnomadttn$START_hg38[v] & ttn_ex$END>=gnomadttn$END_hg38[v]) | (ttn_ex$START>=gnomadttn$START_hg38[v] & ttn_ex$START<=gnomadttn$END_hg38[v] & ttn_ex$END>=gnomadttn$END_hg38[v]) |  (ttn_ex$END<=gnomadttn$END_hg38[v] & ttn_ex$END>=gnomadttn$START_hg38[v] & ttn_ex$START<=gnomadttn$START_hg38[v]) | (ttn_ex$START>=gnomadttn$START_hg38[v] & ttn_ex$END<=gnomadttn$END_hg38[v])))
  if(n_ex>0){
    tokeep = c(tokeep,T)
  }
  else{
    tokeep = c(tokeep,F)
  }
}

gnomadttn = gnomadttn[as.logical(tokeep),]
gnomad = rbind(gnomadnottn,gnomadttn)

gnomad$CHROM = as.numeric(as.character(gsub("chr","",gnomad$CHROM)))
gnomad = gnomad[order(gnomad$CHROM,gnomad$START_hg38,decreasing=F),]

gnomad$CHROM = paste("chr",as.character(gnomad$CHROM),sep="")


#####################################
######process exomedepth results#####
#####################################

#add exon information to the gnomad sv callset
add_exon_info(gnomad,"gnomad")

#remove variants not hitting any exon, and remove inversions
gnomad = subset(gnomad,(gnomad$EXONS != "" & !grepl("inversion",gnomad$Consequence)))

#output gnomAD burden figures
denom = round(mean(gnomad$Allele.Number)/2) #divided by 2 as we count carriers
num = sum(sapply(1:nrow(gnomad), function(x) gnomad$Allele.Count[x]-gnomad$Homozygote.Count[x])) #count carriers
pc = num/denom*100
print(paste("The burden of large coding variants (>50bp) in gnomAD in the analysed genes is ",num,"/",denom," (",round(pc,digits=2),"%). To be conservative and to ensure we don't exclude any TP, we derive a BF cutoff based on a burden of 0.5% for these variants in HVOL. ",sep=""))

#filter HVOL variant calls to keep the same burden of variants, prioritizing variants by higher BF values
n_vars_to_keep = round(902/200) #this keeps the top-BF variants in HVOL to make up for a 0.5% burden
ed_hvol_final = ed_hvol[order(ed_hvol$BF,decreasing = T),]
ed_hvol_discarded = ed_hvol_final[(n_vars_to_keep+1):nrow(ed_hvol_final),]
ed_hvol_final = ed_hvol_final[1:n_vars_to_keep,]
ed_hvol_final = ed_hvol_final[order(ed_hvol_final$chromosome,ed_hvol_final$start.p),]

#derive BF cutoff to be applied to DCM and HCM variants, and apply it
bf_cutoff = min(ed_hvol_final$BF)

print(paste("This corresponds to a BF cutoff of BF>",bf_cutoff,", which is applied to ExomeDepth calls in all cohorts",sep=""))

ed_dcm_final = subset(ed_dcm,ed_dcm$BF>bf_cutoff)
ed_dcm_discarded = subset(ed_dcm,ed_dcm$BF<=bf_cutoff)
ed_hcm_final = subset(ed_hcm,ed_hcm$BF>bf_cutoff)
ed_hcm_discarded = subset(ed_hcm,ed_hcm$BF<=bf_cutoff)

#create a single entry for each recurrent variant in ed results
fix_ed_results_format(ed_dcm_final,"ed_dcm_final")
fix_ed_results_format(ed_hcm_final,"ed_hcm_final")
fix_ed_results_format(ed_hvol_final,"ed_hvol_final")
fix_ed_results_format(ed_dcm_discarded,"ed_dcm_discarded")
fix_ed_results_format(ed_hcm_discarded,"ed_hcm_discarded")
fix_ed_results_format(ed_hvol_discarded,"ed_hvol_discarded")

#save lists of variants discarded by the BF filter
write.table(ed_hvol_discarded,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_HVOL_l_RESULTS_ttnhighpsi_discardedByBF.tsv",row.names=F,col.names=T,sep="\t",quote=F)
write.table(ed_hcm_discarded,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_HCM_l_RESULTS_ttnhighpsi_discardedByBF.tsv",row.names=F,col.names=T,sep="\t",quote=F)
write.table(ed_dcm_discarded,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/exomedepth_DCM_l_RESULTS_ttnhighpsi_discardedByBF.tsv",row.names=F,col.names=T,sep="\t",quote=F)

#remove old variant tables
rm(ed_dcm,ed_hcm,ed_hvol)

#####################################
###process graf/gatk/manta results###
#####################################

#merge the various callers' results formats
#convert POS to START and END in the relevant tables, to enable the addition of exon info
add_start_end_gatk_graf(graf_dcm,"graf_dcm")
add_start_end_gatk_graf(graf_hcm,"graf_hcm")
add_start_end_gatk_graf(graf_hvol_h,"graf_hvol_h")
add_start_end_gatk_graf(graf_hvol_d,"graf_hvol_d")
add_start_end_gatk_graf(gatk_hcm,"gatk_hcm")
add_start_end_gatk_graf(gatk_dcm,"gatk_dcm")
add_start_end_gatk_graf(gatk_hvol,"gatk_hvol")
add_start_end_manta(manta_dcm,"manta_dcm")

#add exon info also to these results and fix column order
add_exon_info(graf_dcm,"graf_dcm")
add_exon_info(graf_hcm,"graf_hcm")
add_exon_info(graf_hvol_h,"graf_hvol_h")
add_exon_info(graf_hvol_d,"graf_hvol_d")
add_exon_info(gatk_dcm,"gatk_dcm")
add_exon_info(gatk_hcm,"gatk_hcm")
add_exon_info(gatk_hvol,"gatk_hvol")
add_exon_info(manta_dcm,"manta_dcm")


####################################
##add caller/cohort info and merge##
####################################
add_column(ed_dcm_final,"ed_dcm_final","COHORT","DCM")
add_column(ed_dcm_final,"ed_dcm_final","CALLER","ED")
add_column(ed_hcm_final,"ed_hcm_final","COHORT","HCM")
add_column(ed_hcm_final,"ed_hcm_final","CALLER","ED")
add_column(ed_hvol_final,"ed_hvol_final","COHORT","HVOL")
add_column(ed_hvol_final,"ed_hvol_final","CALLER","ED")
add_column(graf_dcm,"graf_dcm","COHORT","DCM")
add_column(graf_dcm,"graf_dcm","CALLER","GRAF")
add_column(graf_hcm,"graf_hcm","COHORT","HCM")
add_column(graf_hcm,"graf_hcm","CALLER","GRAF")
add_column(graf_hvol_d,"graf_hvol_d","COHORT","HVOL")
add_column(graf_hvol_d,"graf_hvol_d","CALLER","GRAF")
add_column(graf_hvol_h,"graf_hvol_h","COHORT","HVOL")
add_column(graf_hvol_h,"graf_hvol_h","CALLER","GRAF")
add_column(gatk_dcm,"gatk_dcm","COHORT","DCM")
add_column(gatk_dcm,"gatk_dcm","CALLER","GATK")
add_column(gatk_hcm,"gatk_hcm","COHORT","HCM")
add_column(gatk_hcm,"gatk_hcm","CALLER","GATK")
add_column(gatk_hvol,"gatk_hvol","COHORT","HVOL")
add_column(gatk_hvol,"gatk_hvol","CALLER","GATK")
add_column(manta_dcm,"manta_dcm","COHORT","DCM")
add_column(manta_dcm,"manta_dcm","CALLER","MANTA")

#FROM HERE THE VARIANT FILTERING IS MANUAL, AS THERE ARE TOO MANY SPECIFIC CASES TO BE HANDLED AUTOMATICALLY#
#graf_hvol_d and graf_hvol_h contain exactly the same variants, with the only difference being the n_samples 
#in which a recurrent ACTC1 variant is called - from now on just proceed with a single graf_hvol variant table
graf_hvol = graf_hvol_d
rm(graf_hvol_d,graf_hvol_h)

#concatenate all tables into one, to be analyzed manually
final_candidates = rbind(graf_hvol,graf_dcm,graf_hcm,gatk_hvol,gatk_dcm,gatk_hcm,manta_dcm,ed_dcm_final,ed_hcm_final,ed_hvol_final)

#remove the chr prefix (gatk/graf/manta), adapt type nomenclature, fix SVLEN and fix n_samples to show the n of carriers
final_candidates$CHROM = gsub("chr","",final_candidates$CHROM)
final_candidates$TYPE = gsub("deletion","loss",final_candidates$TYPE)
final_candidates$TYPE = gsub("duplication","gain",final_candidates$TYPE)
final_candidates$SVLEN = as.numeric(as.character(final_candidates$SVLEN))
final_candidates$SVLEN = abs(final_candidates$SVLEN)
final_candidates$N_SAMPLES = sapply(1:nrow(final_candidates),function(x) length(unlist(strsplit(final_candidates$SAMPLES[x],","))))
final_candidates = final_candidates[order(as.numeric(final_candidates$CHROM),as.numeric(final_candidates$START_hg38),as.numeric(final_candidates$END_hg38)),]
final_candidates[which(is.na(final_candidates$gnomAD_AF)),"gnomAD_AF"] = 0
VAR_ID = 1:nrow(final_candidates)
final_candidates = cbind(VAR_ID,final_candidates)

#save table
write.table(final_candidates,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/FINAL_CANDIDATES.tsv",row.names=F,col.names=T,sep="\t",quote=F)
