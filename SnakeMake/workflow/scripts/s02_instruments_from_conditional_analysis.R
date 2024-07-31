suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of Cojo"),
  make_option("--conditional_output", default=NULL, help="Output path and name for list of instruments from Cojo"),
  make_option("--unconditional_output", default=NULL, help="Output path and name for list of instruments from Cojo unconditional")
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cojo<-fread(opt$input)
cond_path<-opt$cond_output
uncond_path<-opt$uncond_output
mapping<-fread(opt$mapping)

##mapping file
# mapping<-fread("/home/solene.cadiou/basic_GWAS_protein/meta_results/MR/MR_instruments_selection/mapped_gene_file_GRCh37_21052025.txt")
mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)
##################
#add print before and after for

###cis mapping of Cojo file
for (i in 1:nrow(cojo)){
  cojo$cis_or_trans[i]<-map(cojo$study_id[i],cojo$Chr[i],cojo$bp[i],cojo$bp[i],mapping)
}
##filtering only cis
cojo<-cojo[which(cojo$cis_or_trans=="cis"),]

##Filter p-value bonferroni corrected
cojo_sign<-cojo[which(cojo$p < (5e-8/3958)), ]
cojo_j_sign<- cojo[which(cojo$pJ < (5e-8/3958)), ]
##Fstats computation and filtering
cojo_sign$Fstats<-((cojo_sign$b^2)/(cojo_sign$se^2))
cojo_sign<-cojo_sign[which(cojo_sign$Fstats>=10), ]

cojo_j_sign$Fstats_j<-((cojo_j_sign$bJ^2)/(cojo_j_sign$bJ_se^2))
cojo_j_sign<-cojo_j_sign[which(cojo_j_sign$Fstats_j>=10), ]

write.table(cojo_sign,uncond_path)
write.table(cojo_j_sign, cond_path)


