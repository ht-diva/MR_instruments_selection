suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of Cojo"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--conditional_output", default=NULL, help="Output path and name for list of instruments from Cojo"),
  make_option("--unconditional_output", default=NULL, help="Output path and name for list of instruments from Cojo unconditional"),
  make_option("--cojo_output", default=NULL, help="Output path from Cojo with cis and trans"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cojo<-fread(opt$input)
conditional_path<-opt$conditional_output
unconditional_path<-opt$unconditional_output
cojo_path<-opt$cojo_output
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
cojo_cis<-cojo[which(cojo$cis_or_trans=="cis"),]
##select unconditional beta se
cojo_cis$Fstats_C<-((cojo_cis$bC^2)/(cojo_cis$bC_se^2))
cojo_cis<-cojo_cis[which(cojo_cis$Fstats_C>=10), ]

#From mapping file we have to save the following header:
#GENE_NAMEDATASET >
#TISSUE
#SNPID
#RSID
#CHR > Chr
#POS_38
#POS_37 BETA
#SE
#MinusLog10PVAL
#EFFECT_ALLELE
#OTHER_ALLELE
#MAF
#EAF
#SAMPLESIZE
#GENE_ENSEMBL
#TSS_37
#TSS_38
#SeqID
#UNIPROT
#PROTEIN_NAME
#PROTEIN_LONG_NAME
#PVE
#FILENAME
#Gene.type

#we have to save the documents as MVP wants:

cojo_conditional<-cojo_cis[, c("Chr","SNP","bp","refA","freq","n","freq_geno",
                        "bC","bC_se","pC","LD_r","snp_map","sdY","study_id","locus","EA", "Fstats_j")]

cojo_unconditional<-cojo_cis[, c("Chr","SNP","bp","refA","freq","b","se","p","n","freq_geno",
                        "LD_r","snp_map","sdY","study_id","locus","EA","Fstats_j")]

fwrite(cojo_conditional, conditional_path)
fwrite(cojo_unconditional,unconditional_path)
fwrite(cojo, cojo_path)