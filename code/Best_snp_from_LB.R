#install.packages(data.table)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("IRanges")
library(data.table)
library(IRanges)
source("/home/solene.cadiou/basic_GWAS_protein/meta_results/MR/MR_instruments_selection/code/mapping_function_for_regions.R")

###loading and parameters########
##Locus_breaker results
LB<-fread("/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats/24-May-24_regional_associations_with_6236_proteins_mhc_excluded.csv")
##number of effective tests
NEF<-3978
##path to sumstats
path_to_sumstats<-"/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats_digits/output/"

##mapping file
mapping<-fread("/home/solene.cadiou/basic_GWAS_protein/meta_results/MR/MR_instruments_selection/code/mapped_gene_file_GRCh37_21052025.txt")
mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)
##################

###cis mapping of LB file
# LB$cis_or_trans<-apply(LB[,c("study_id","chr","start","end")],1,function(X) map(seqId=X["study_id"],chr=X["chr"],start=as.numeric(X["start"]),end=as.numeric(X["end"]),mapping_file=mapping))
for (i in 1:nrow(LB)){
    LB$cis_or_trans[i]<-map(LB$study_id[i],LB$chr[i],LB$start[i],LB$end[i],mapping)
  } ##debug the apply command instead of using loop
table(LB$cis_or_trans)

##filtering only cis
LB<-LB[LB$cis_or_trans=="cis",]

##Fstats computation
LB$Fstats<-((LB$BETA^2)/(LB$SE^2))
##filtering snp passing fstats threshold
LB$instrum<-LB$Fstats>=10 ##to check with Giulia >10 or >=10
# table(LB$instrum)

##selecting only columns of interests
LB<-LB[,c("start","end","chr","POS","SNPID","EA","NEA","EAF","BETA","SE","MLOG10P","Fstats","instrum","study_id","cis_or_trans")]
##here all the rows with instrum =T are instruments
##TO CHECK WITH MVP: if 2 rows with identical seqID, what should we do? keep both? choose the one with strongest pval?
############################################

##we need to define instruments from the sunmstats for the regions in which the top snp does not pass Fstats filter
if (FALSE%in%LB$instrum){
index_not_passing<-which(!LB$instrum)
LB_not_passing<-LB[!LB$instrum,]
for (i in 1:nrow(LB_not_passing)){
  ##build path (to modify with the right path)
  path<-paste("/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats_digits/output/",LB_not_passing$study_id[i],"/",LB_not_passing$study_id[i],".gwaslab.tsv.gz",sep="")
  sumstats<-fread(path,select = c("CHR","POS","SNPID","EA","NEA","EAF","BETA","SE","MLOG10P"))
  cis<-sumstats[sumstats$CHR==LB_not_passing$chr[i],]
  cis<-cis[cis$POS>=LB_not_passing$start[i]&cis$POS<=LB_not_passing$end[i],]
  cis<-cis[cis$MLOG10P>=(-log10(5/NEF*10^(-8))),]
  cis$Fstats<-(cis$BETA^2)/(cis$SE^2)
  cis<-cis[Fstats>=10,]
  if (nrow(cis)>=1){
    cis<-cis[which.max(cis$MLOG10P)]
    cis$instrum<-TRUE
    ##check
    table(colnames(cis)[-1]==colnames(LB_not_passing)[4:13])
    ##check
    LB_not_passing$study_id[i]==LB$study_id[index_not_passing[i]]
    LB[index_not_passing[i],]<-c(NA,NA,cis[1,], LB_not_passing$study_id[i],"cis")
  }
}
}


##selecting only columns of interests
LB<-LB[,c("chr","POS","SNPID","EA","NEA","EAF","BETA","SE","MLOG10P","study_id","cis_or_trans","Fstats","instrum")]
##homogenization with meta-analysis colnames
colnames(LB)[1]<-"CHR"
##save
write.table(LB,"best_snps_from_LB.txt") ##add path


##########################
##TESTS##################
map(seqId="seq.10025.1",chr=LB$chr[LB$study_id=="seq.10025.1"],start=as.numeric(LB$start[LB$study_id=="seq.10025.1"]),
    end=as.numeric(LB$end[LB$study_id=="seq.10025.1"]),mapping_file=mapping)

LB[3,]$instrum<-FALSE
LB[5,]$instrum<-FALSE
#######################
