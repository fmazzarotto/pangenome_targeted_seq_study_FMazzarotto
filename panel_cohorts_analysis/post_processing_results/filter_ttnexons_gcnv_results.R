filter_ttn_exons_gcnv <- function(samplefile,cohort){
  genes = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/genelists/genes.tsv",header=T)
  
  id = gsub(".tsv","",basename(samplefile))
  vars = read.table(samplefile,header=T,sep="\t")[,c(1:2,4:5,8)]
  ttnhpexons = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/ttnexons/TTN_highPSI_GRCh38.bed",header=F)
  gnomad = read.csv("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/gnomad/gnomAD_SV_CMgenes.csv",header=T,colClasses = c(rep("numeric",5),rep("character",3),rep("numeric",6)))
  
  ttnvars = subset(vars,vars$CHROM == "chr2" & vars$POS >= 178525989 & vars$POS <= 178807423)
  nottnvars = subset(vars,!(vars$CHROM == "chr2" & vars$POS >= 178525989 & vars$POS <= 178807423))
  
  if(sum(nrow(ttnvars),nrow(nottnvars)) != nrow(vars)){
    stop("There is something wrong in how TTN and other genes' variants are extracted - check code")    
  }
  
  tokeep = logical(nrow(ttnvars))
  
  if(nrow(ttnvars)>0){
    for(v in 1:nrow(ttnvars)){
      #extract variant start and end coordinates
      v_start = as.numeric(as.character(unlist(strsplit(ttnvars$ID[v],"_"))[3]))
      v_end = as.numeric(as.character(unlist(strsplit(ttnvars$ID[v],"_"))[4]))
    
      #create a vector with overlap T or F for each high-psi exon relative to the variant
      overlap = (ttnhpexons[,2] <= v_end & ttnhpexons[,3] >= v_start)
      
      #if any exon has overlap=T. store T in the variant's element of tokeep
      tokeep[v] = any(overlap)
    }
  
    ttnvars = ttnvars[as.logical(tokeep),]
    newvars = rbind(ttnvars,nottnvars)
  } else{
    newvars = nottnvars
  }
  
  if(nrow(newvars) > 0){
    #extract QUAL (not doable with tableize.vcf.py and integrate it into the tsv)
    vcfvars = read.table(text = system(paste0("zless /rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/",cohort,"/",id,".FINAL.vcf.gz | grep -v ^#"),intern=T), sep="\t", header=F, stringsAsFactors=F)[,c(1:3,5:6)]
    colnames(vcfvars) = c("CHROM","POS","ID","ALT","QUAL")
    newvars = merge(newvars,vcfvars,by=c("CHROM","POS","ID","ALT"),all.x=T)
    newvars = unique(newvars)
    #extract start and end pos of each variant
    newvars$START_hg38 = as.numeric(sapply(1:nrow(newvars),function(x) unlist(strsplit(as.character(newvars$ID[x]),"_"))[3]))
    newvars$END_hg38 = as.numeric(sapply(1:nrow(newvars),function(x) unlist(strsplit(as.character(newvars$ID[x]),"_"))[4]))
    newvars$SVLEN = newvars$END_hg38 - newvars$START_hg38
    
    
    #add exon info to each variant
    library(biomaRt)
    ensmart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    target_bed = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/TSC_40bpOverhang_hg38.bed")
    colnames(target_bed) = c("chromosome","start","end","name")
    target_bed$name = gsub("SEPN1","SELENON",target_bed$name) #needed to avoid NAs in fetching strand info
    target_bed$name = gsub("MURC","CAVIN4",target_bed$name)
    target_bed$name = gsub("TAZ","TAFAZZIN",target_bed$name)
    genestrandinfo = getBM(attributes=c("hgnc_symbol","strand"),filters="hgnc_symbol",values=target_bed$name,mart=ensmart)
    target_bed$strand = genestrandinfo$strand[match(target_bed$name, genestrandinfo$hgnc_symbol)]
    exon_id = c()
    for(g in unique(target_bed$name)){
      ei = c()
      gs = subset(target_bed,target_bed$name==g)
      if(gs$strand[1] == 1){
        for(i in 1:nrow(gs)){
          ei = c(ei,paste(g,"_",i,sep=""))
        }
      }
      else{
        exon_no = nrow(gs)
        for(i in 1:nrow(gs)){
          ei = c(ei,paste(g,"_",exon_no-i+1,sep=""))
        }
      }
      exon_id = c(exon_id,ei)
    }
    target_bed = cbind(target_bed,exon_id)
    EXONS = c()
    GENES = c()
    EXON_RANGE = c()
    ENSG_ID = c()
    ENST_ID = c()
    for(cnv in 1:nrow(newvars)){
      exonsinvolved = subset(target_bed,(target_bed$chromosome==newvars$CHROM[cnv] & ((target_bed$start<=newvars$START_hg38[cnv] & target_bed$end>=newvars$START_hg38[cnv]) | (target_bed$start>=newvars$START_hg38[cnv] & target_bed$end<=newvars$END_hg38[cnv]) | (target_bed$start<=newvars$START_hg38[cnv] & target_bed$end>=newvars$END_hg38[cnv]) | (target_bed$start<=newvars$START_hg38[cnv] & target_bed$end>=newvars$END_hg38[cnv]))))
      exonsinvolved = subset(exonsinvolved,exonsinvolved$name %in% genes$GENE)
      exonids = paste(exonsinvolved$exon_id,collapse=",")
      EXONS = c(EXONS,exonids)
      GENES = c(GENES,paste(unique(exonsinvolved$name),collapse=","))
      
      if(length(unique(exonsinvolved$name))==1){
        exon_numbers = as.numeric(sub(".*_", "", exonsinvolved$exon_id))
        exonrange = paste0(min(exon_numbers), "-", max(exon_numbers))
        EXON_RANGE = c(EXON_RANGE,exonrange)
        ENSG_ID = c(ENSG_ID,genes[which(genes$GENE == unique(exonsinvolved$name)),"ENSG"])
        ENST_ID = c(ENST_ID,genes[which(genes$GENE == unique(exonsinvolved$name)),"ENST"])
      }
      else{
        EXON_RANGE = c(EXON_RANGE,"Check manually")
        ENSG_ID = c(ENSG_ID,"Check manually")
        ENST_ID = c(ENST_ID,"Check manually")
      }
    }
    newvars = cbind(newvars,GENES,EXONS,EXON_RANGE,ENSG_ID,ENST_ID)
    
    #add exon info to gnomAD table too
    EXONS = c()
    gnomad$CHROM = paste0("chr",gnomad$CHROM)
    for(cnv in 1:nrow(gnomad)){
      exonsinvolved = subset(target_bed,(target_bed$chromosome==gnomad$CHROM[cnv] & ((target_bed$start<=gnomad$START_hg38[cnv] & target_bed$end>=gnomad$START_hg38[cnv]) | (target_bed$start>=gnomad$START_hg38[cnv] & target_bed$end<=gnomad$END_hg38[cnv]) | (target_bed$start<=gnomad$START_hg38[cnv] & target_bed$end>=gnomad$END_hg38[cnv]) | (target_bed$start<=gnomad$START_hg38[cnv] & target_bed$end>=gnomad$END_hg38[cnv]))))
      exonids = paste(exonsinvolved$exon_id,collapse=",")
      EXONS = c(EXONS,exonids)
    }
    gnomad = cbind(gnomad,EXONS)
    
    gnomad_common = subset(gnomad,gnomad$Allele.Frequency >= 0.0001)
    
    #add/check gnomAD AF
    gnomAD_AF = c()
    for(v in 1:nrow(newvars)){
      
      if(sum(gnomad$EXONS == newvars$EXONS[v]) > 0){
        gnomAD_AF = c(gnomAD_AF,"Check_manually")
      }
      else{
        gnomAD_AF = c(gnomAD_AF,0)
      }
    }
    newvars = cbind(newvars,gnomAD_AF)
    
    
    REF = rep(NA,nrow(newvars))
    REF_orig = rep("N",nrow(newvars))
    ALT = rep(NA,nrow(newvars))
    ED_BF = rep(NA,nrow(newvars))
    CSQ = rep(NA,nrow(newvars))
    HGVSc = rep(NA,nrow(newvars))
    HGVSp = rep(NA,nrow(newvars))
    CALLER = rep("GATKgCNV",nrow(newvars))
    TYPE = ifelse(newvars$ALT == "<DEL>","loss","gain")
    COHORT = rep(cohort,nrow(newvars))
    
    colnames(newvars)[2] = "POS_orig"
    colnames(newvars)[4] = "ALT_orig"
    
    
    newvars = cbind(newvars,REF,REF_orig,ALT,ED_BF,CSQ,HGVSc,HGVSp,CALLER,TYPE,ENSG_ID,ENST_ID,gnomAD_AF)
    newvars = newvars[,c("CHROM","START_hg38","END_hg38","REF","ALT","TYPE","QUAL","SVLEN","ED_BF","POS_orig","REF_orig","ALT_orig","SAMPLES","GENES","ENSG_ID","ENST_ID","CSQ","HGVSc","HGVSp","gnomAD_AF","EXON_RANGE")]

    write.table(newvars,paste0("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/",cohort,"/",id,"_TTNhp.tsv"),col.names=T,row.names=F,sep="\t",quote=F)
  }
  else{
    dir.create(paste0("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/",cohort,"/ttn_lowpsi"),showWarnings = F)
    file.rename(paste0("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/",cohort,"/",id,".FINAL.vcf.gz"),paste0("/rds/general/project/lms-ware-analysis/live/francesco/pangenomes_cnv/gatk_gcnv_analysis/results/",cohort,"/ttn_lowpsi/",id,".FINAL.vcf.gz"))
  }
  file.remove(samplefile)
}

#main
args = commandArgs(trailingOnly = TRUE)
sample_file <- args[1]
cohort <- args[2]
filter_ttn_exons_gcnv(sample_file,cohort)
