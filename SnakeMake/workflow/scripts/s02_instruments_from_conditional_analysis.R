suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of Cojo"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--liftover", default=NULL, help="liftover file path"),
  make_option("--conditional_output", default=NULL, help="Output path and name for list of instruments from Cojo"),
  make_option("--unconditional_output", default=NULL, help="Output path and name for list of instruments from Cojo unconditional"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cojo<-fread(opt$input)
mapping<-fread(opt$mapping)
liftover<-fread(liftover$mapping)
conditional_path<-opt$conditional_output
unconditional_path<-opt$unconditional_output

cojo <- cojo[order(cojo$Chr, cojo$bp), ]
split_variants <- strsplit(cojo$SNP, ":")
cojo$EA <- sapply(split_variants, function(x) x[3])
cojo$NEA <- sapply(split_variants, function(x) x[4])

liftover <- liftover[,c(1:3)]
cojo <- cbind(liftover, cojo)

mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)

for (i in 1:nrow(cojo)){
  cojo$cis_or_trans[i]<-map(cojo$study_id[i],cojo$Chr[i],cojo$bp[i],cojo$bp[i],mapping)
}

cojo_cis<-cojo[which(cojo$cis_or_trans=="cis"),]
##select unconditional beta se
cojo_cis$Fstats_conditional<-((cojo_cis$bC^2)/(cojo_cis$bC_se^2))
cojo_cis_conditional<-cojo_cis[which(cojo_cis$Fstats_conditional>=10), ]
cojo_cis$Fstats<-((cojo_cis$b^2)/(cojo_cis$se^2))
cojo_cis_unconditional <- cojo_cis[which(cojo_cis$Fstats>=10), ]

merged_unconditional <- cojo_cis_unconditional %>%
  left_join(mapping, by = c("study_id" = "target"), relationship = "many-to-many") %>%
  filter(Chr == chromosome, (bp >= cis_start & bp <= cis_end))

merged_unconditional$DATASET="INTERVAL_CHRIS_META_COJO"
merged_unconditional$TISSUE="WholeBlood"
merged_unconditional$FILENAME=NA
merged_unconditional$Gene.type = "protein_coding"
merged_unconditional$MAF <- pmin(merged_unconditional$freq_geno, 1-merged_unconditional$freq_geno)
merged_unconditional$PVE <- (2*(merged_unconditional$b^2)*merged_unconditional$MAF*(1-merged_unconditional$MAF))/
  (2*(merged_unconditional$b^2)*merged_unconditional$MAF*(1-merged_unconditional$MAF)+(merged_unconditional$se^2)*2*merged_unconditional$n*merged_unconditional$MAF*(1-merged_unconditional$MAF))

merged_unconditional <- merged_unconditional %>%
  dplyr::select(DATASET, TISSUE, SNP, Chr, bp, V2, b, se, mlog10p, EA, NEA, MAF,
                freq_geno, n, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, study_id, UniProt_ID,
                Target_Name, Target_Full_Name, PVE, FILENAME, Gene.type)

names(merged_unconditional)[names(merged_unconditional) == "b"] <- "BETA"
names(merged_unconditional)[names(merged_unconditional) == "se"] <- "SE"
names(merged_unconditional)[names(merged_unconditional) == "bp"] <- "POS_37"
names(merged_unconditional)[names(merged_unconditional) == "V2"] <- "POS_38"
names(merged_unconditional)[names(merged_unconditional) == "freq_geno"] <- "EAF"
names(merged_unconditional)[names(merged_unconditional) == "MLOG10P"] <- "MinusLog10PVAL"
names(merged_unconditional)[names(merged_unconditional) == "EA"] <- "EFFECT_ALLELE"
names(merged_unconditional)[names(merged_unconditional) == "NEA"] <- "OTHER_ALLELE"
names(merged_unconditional)[names(merged_unconditional) == "n"] <- "SAMPLESIZE"
names(merged_unconditional)[names(merged_unconditional) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged_unconditional)[names(merged_unconditional) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged_unconditional)[names(merged_unconditional) == "TSS"] <- "TSS_37"
names(merged_unconditional)[names(merged_unconditional) == "study_id"] <- "SeqID"
names(merged_unconditional)[names(merged_unconditional) == "UniProt_ID"] <- "UNIPROT"
names(merged_unconditional)[names(merged_unconditional) == "Target_Name"] <- "PROTEIN_NAME"
names(merged_unconditional)[names(merged_unconditional) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

cojo_unconditional<-merged_unconditional

merged_conditional <- cojo_cis_conditional %>%
  left_join(mapping, by = c("study_id" = "target"), relationship = "many-to-many") %>%
  filter(Chr == chromosome, (bp >= cis_start & bp <= cis_end))

merged_conditional$DATASET="INTERVAL_CHRIS_META_COJO"
merged_conditional$TISSUE="WholeBlood"
merged_conditional$FILENAME=NA
merged_conditional$Gene.type = "protein_coding"
merged_conditional$MAF <- pmin(merged_conditional$freq_geno, 1-merged_conditional$freq_geno)
merged_conditional$PVE_conditional <- (2*(merged_conditional$bC^2)*merged_conditional$MAF*(1-merged_conditional$MAF))/
  (2*(merged_conditional$bC^2)*merged_conditional$MAF*(1-merged_conditional$MAF)+(merged_conditional$bC_se^2)*2*merged_conditional$n*merged_conditional$MAF*(1-merged_conditional$MAF))

merged_conditional <- merged_conditional %>%
  dplyr::select(DATASET, TISSUE, SNP, Chr, bp, V2, bC, bC_se, mlog10pC, EA, NEA, MAF,
                freq_geno, n, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, study_id, UniProt_ID,
                Target_Name, Target_Full_Name, PVE_conditional, FILENAME, Gene.type)

names(merged_conditional)[names(merged_conditional) == "bC"] <- "BETA"
names(merged_conditional)[names(merged_conditional) == "bC_se"] <- "SE"
names(merged_conditional)[names(merged_conditional) == "bp"] <- "POS_37"
names(merged_conditional)[names(merged_conditional) == "V2"] <- "POS_38"
names(merged_conditional)[names(merged_conditional) == "freq_geno"] <- "EAF"
names(merged_conditional)[names(merged_conditional) == "PVE_conditional"] <- "PVE"
names(merged_conditional)[names(merged_conditional) == "mlog10pC"] <- "MinusLog10PVAL"
names(merged_conditional)[names(merged_conditional) == "EA"] <- "EFFECT_ALLELE"
names(merged_conditional)[names(merged_conditional) == "NEA"] <- "OTHER_ALLELE"
names(merged_conditional)[names(merged_conditional) == "n"] <- "SAMPLESIZE"
names(merged_conditional)[names(merged_conditional) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged_conditional)[names(merged_conditional) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged_conditional)[names(merged_conditional) == "TSS"] <- "TSS_37"
names(merged_conditional)[names(merged_conditional) == "study_id"] <- "SeqID"
names(merged_conditional)[names(merged_conditional) == "UniProt_ID"] <- "UNIPROT"
names(merged_conditional)[names(merged_conditional) == "Target_Name"] <- "PROTEIN_NAME"
names(merged_conditional)[names(merged_conditional) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

cojo_conditional<-merged_conditional

fwrite(cojo_conditional, conditional_path)
fwrite(cojo_unconditional,unconditional_path)

