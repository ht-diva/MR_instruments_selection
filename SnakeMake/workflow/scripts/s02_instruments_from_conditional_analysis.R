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
liftover<-fread(opt$liftover)
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

collapsed_df_unconditional <- merged_unconditional %>%
  group_by(study_id, SNP) %>%
  # Step 3: Collapse all columns by concatenating values where necessary
  summarise(
    across(everything(), ~ paste(unique(.), collapse = "|"), .names = "collapsed_{col}"),
    .groups = "drop")


collapsed_df_unconditional$DATASET="INTERVAL_CHRIS_META_COJO"
collapsed_df_unconditional$TISSUE="WholeBlood"
collapsed_df_unconditional$FILENAME=NA
collapsed_df_unconditional$Gene.type = "protein_coding"
collapsed_df_unconditional$N <- 13445
collapsed_df_unconditional$collapsed_freq_geno <- as.numeric(collapsed_df_unconditional$collapsed_freq_geno)
collapsed_df_unconditional$collapsed_b <- as.numeric(collapsed_df_unconditional$collapsed_b)
collapsed_df_unconditional$collapsed_se <- as.numeric(collapsed_df_unconditional$collapsed_se)

collapsed_df_unconditional$MAF <- pmin(collapsed_df_unconditional$collapsed_freq_geno, 1-collapsed_df_unconditional$collapsed_freq_geno)
collapsed_df_unconditional$PVE <- (2*(collapsed_df_unconditional$collapsed_b^2)*collapsed_df_unconditional$MAF*(1-collapsed_df_unconditional$MAF))/
  (2*(collapsed_df_unconditional$collapsed_b^2)*collapsed_df_unconditional$MAF*(1-collapsed_df_unconditional$MAF)+(collapsed_df_unconditional$collapsed_se^2)*2*collapsed_df_unconditional$N*collapsed_df_unconditional$MAF*(1-collapsed_df_unconditional$MAF))

collapsed_df_unconditional <- collapsed_df_unconditional %>%
  dplyr::select(DATASET, TISSUE, SNP, collapsed_Chr, collapsed_bp, collapsed_V2, collapsed_locus_START_END_37,
                collapsed_b, collapsed_se, collapsed_mlog10p, collapsed_EA, collapsed_NEA, MAF,
                collapsed_freq_geno, N, PVE,collapsed_Entrez_Gene_Name, collapsed_Ensembl_Gene_ID,
                collapsed_TSS, study_id, collapsed_UniProt_ID,
                collapsed_Target_Name, collapsed_Target_Full_Name, FILENAME, Gene.type)


names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Chr"] <- "CHR"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_bp"] <- "POS_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_V2"] <- "POS_38"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_b"] <- "BETA"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_se"] <- "SE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_mlog10p"] <- "MinusLog10PVAL"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_EA"] <- "EFFECT_ALLELE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_NEA"] <- "OTHER_ALLELE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_freq_geno"] <- "EAF"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_N"] <- "SAMPLESIZE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Entrez_Gene_Name"] <- "GENE_NAME"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_TSS"] <- "TSS_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_study_id"] <- "SeqID"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_UniProt_ID"] <- "UNIPROT"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Target_Name"] <- "PROTEIN_NAME"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Target_Full_Name"] <- "PROTEIN_LONG_NAME"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Gene.type"] <- "Gene.type"

cojo_unconditional<-collapsed_df_unconditional

merged_conditional <- cojo_cis_conditional %>%
  left_join(mapping, by = c("study_id" = "target"), relationship = "many-to-many") %>%
  filter(Chr == chromosome, (bp >= cis_start & bp <= cis_end))

collapsed_df_conditional <- merged_conditional %>%
  group_by(study_id, SNP) %>%
  # Step 3: Collapse all columns by concatenating values where necessary
  summarise(
    across(everything(), ~ paste(unique(.), collapse = "|"), .names = "collapsed_{col}"),
    .groups = "drop")


collapsed_df_conditional$DATASET="INTERVAL_CHRIS_META_COJO"
collapsed_df_conditional$TISSUE="WholeBlood"
collapsed_df_conditional$FILENAME=NA
collapsed_df_conditional$Gene.type = "protein_coding"
collapsed_df_conditional$N <- 13445
collapsed_df_conditional$collapsed_freq_geno <- as.numeric(collapsed_df_conditional$collapsed_freq_geno)
collapsed_df_conditional$collapsed_bC <- as.numeric(collapsed_df_conditional$collapsed_bC)
collapsed_df_conditional$collapsed_bC_se <- as.numeric(collapsed_df_conditional$collapsed_bC_se)

collapsed_df_conditional$MAF <- pmin(collapsed_df_conditional$collapsed_freq_geno, 1-collapsed_df_conditional$collapsed_freq_geno)
collapsed_df_conditional$PVE <- (2*(collapsed_df_conditional$collapsed_bC_se^2)*collapsed_df_conditional$MAF*(1-collapsed_df_conditional$MAF))/
  (2*(collapsed_df_conditional$collapsed_bC^2)*collapsed_df_conditional$MAF*(1-collapsed_df_conditional$MAF)+(collapsed_df_conditional$collapsed_bC_se^2)*2*collapsed_df_conditional$N*collapsed_df_conditional$MAF*(1-collapsed_df_conditional$MAF))

collapsed_df_conditional <- collapsed_df_conditional %>%
  dplyr::select(DATASET, TISSUE, SNP, collapsed_Chr, collapsed_bp, collapsed_V2, collapsed_locus_START_END_37,
                collapsed_bC, collapsed_bC_se, collapsed_mlog10pC, collapsed_EA, collapsed_NEA, MAF,
                collapsed_freq_geno, N, PVE,collapsed_Entrez_Gene_Name, collapsed_Ensembl_Gene_ID,
                collapsed_TSS, study_id, collapsed_UniProt_ID,
                collapsed_Target_Name, collapsed_Target_Full_Name, FILENAME, Gene.type)


names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Chr"] <- "CHR"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_bp"] <- "POS_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_V2"] <- "POS_38"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_bC"] <- "BETA"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_bC_se"] <- "SE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_mlog10pC"] <- "MinusLog10PVAL"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_EA"] <- "EFFECT_ALLELE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_NEA"] <- "OTHER_ALLELE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_freq_geno"] <- "EAF"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_N"] <- "SAMPLESIZE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Entrez_Gene_Name"] <- "GENE_NAME"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_TSS"] <- "TSS_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_study_id"] <- "SeqID"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_UniProt_ID"] <- "UNIPROT"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Target_Name"] <- "PROTEIN_NAME"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Target_Full_Name"] <- "PROTEIN_LONG_NAME"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Gene.type"] <- "Gene.type"

cojo_conditional<-collapsed_df_conditional

fwrite(cojo_conditional, conditional_path)
fwrite(cojo_unconditional,unconditional_path)