suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
source("mapping_function_for_regions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of Cojo"),
  make_option("--output", default=NULL, help="Output path and name for list of instruments from Cojo"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cojo<-fread(opt$input)
output_path<-opt$output
mapping<-fread(opt$mapping)
####################################
###loading and parameters########
# ##path to sumstats
# path_to_sumstats<-"/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats_digits/output/"

##mapping file
# mapping<-fread("/home/solene.cadiou/basic_GWAS_protein/meta_results/MR/MR_instruments_selection/mapped_gene_file_GRCh37_21052025.txt")
mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)
##################

###cis mapping of Cojo file
# LB$cis_or_trans<-apply(LB[,c("study_id","chr","start","end")],1,function(X) map(seqId=X["study_id"],chr=X["chr"],start=as.numeric(X["start"]),end=as.numeric(X["end"]),mapping_file=mapping))
for (i in 1:nrow(cojo)){
  cojo$cis_or_trans[i]<-map(cojo$study_id[i],cojo$Chr[i],cojo$bp[i],cojo$bp[i],mapping)
}
##filtering only cis
cojo<-cojo[cojo$cis_or_trans=="cis",]

##Fstats computation
cojo$Fstats<-((cojo$b^2)/(cojo$se^2))
cojo$Fstats_j<-((cojo$bJ^2)/(cojo$bJ_se^2))

write.table(cojo,output_path)
