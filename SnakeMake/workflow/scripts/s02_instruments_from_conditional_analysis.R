suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of Cojo"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--input_liftover", default=NULL, help="liftover file path"),
  make_option("--conditional_output", default=NULL, help="Output path and name for list of instruments from Cojo"),
  make_option("--unconditional_output", default=NULL, help="Output path and name for list of instruments from Cojo unconditional"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cojo<-fread(opt$input)
mapping<-fread(opt$mapping)

liftover <- read.table(opt$input_liftover, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
conditional_path<-opt$conditional_output
unconditional_path<-opt$unconditional_output

cojo <- cojo[order(cojo$Chr, cojo$bp), ]
split_variants <- strsplit(cojo$SNP, ":")
cojo$EA <- sapply(split_variants, function(x) x[3])
cojo$NEA <- sapply(split_variants, function(x) x[4])

cojo$locus_START_END_37 <- gsub("^chr[0-9XY]+_", "", cojo$locus)  # Remove "chr" part and chromosome number
cojo$locus_START_END_37 <- gsub("_", "-", cojo$locus_START_END_37)

locus_split <- strsplit(cojo$locus_START_END_37, "-")
cojo$locus_split_start <- sapply(locus_split, function(x) x[1])
cojo$locus_split_end <- sapply(locus_split, function(x) x[2])
cojo$locus_start_extended <- as.numeric(as.character(cojo$locus_split_start))-100000
cojo$locus_start_extended <- ifelse(cojo$locus_start_extende < 0, 0, cojo$locus_start_extende)
cojo$locus_end_extended <- as.numeric(as.character(cojo$locus_split_end))+100000
cojo$locus_extended_START_END_37 <- paste0(cojo$locus_start_extended, "-", cojo$locus_end_extended)

liftover <- liftover[,c(1:3)]
cojo <- cbind(liftover, cojo)

mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)

for (i in 1:nrow(cojo)){
  cojo$cis_or_trans[i]<-map(cojo$study_id[i],cojo$Chr[i],cojo$bp[i],cojo$bp[i],mapping)
}

cojo_cis<-cojo[which(cojo$cis_or_trans=="cis"),]
cojo_cis$N <- 13445
cojo_cis$MAF <- pmin(cojo_cis$freq_geno, 1-cojo_cis$freq_geno)
cojo_cis <- cojo_cis %>%
  group_by(locus_START_END_37) %>%
  mutate(k = n())

cojo_cis$PVE <- (2*(cojo_cis$b^2)*cojo_cis$MAF*(1-cojo_cis$MAF))/
  (2*(cojo_cis$b^2)*cojo_cis$MAF*(1-cojo_cis$MAF)+(cojo_cis$se^2)*2*cojo_cis$N*cojo_cis$MAF*(1-cojo_cis$MAF))
cojo_cis$Fstats <- (cojo_cis$PVE*(cojo_cis$N-2))/(1-cojo_cis$PVE)
cojo_cis$Fstats_multipleMR <- (cojo_cis$PVE*(cojo_cis$N-1-cojo_cis$k))/((1-cojo_cis$PVE)*cojo_cis$k)

cojo_cis$PVE_C <- (2*(cojo_cis$bC^2)*cojo_cis$MAF*(1-cojo_cis$MAF))/
  (2*(cojo_cis$bC^2)*cojo_cis$MAF*(1-cojo_cis$MAF)+(cojo_cis$bC_se^2)*2*cojo_cis$N*cojo_cis$MAF*(1-cojo_cis$MAF))
cojo_cis$Fstats_C <- (cojo_cis$PVE_C*(cojo_cis$N-2))/(1-cojo_cis$PVE_C)
cojo_cis$Fstats_C_multipleMR <- (cojo_cis$PVE_C*(cojo_cis$N-1-cojo_cis$k))/((1-cojo_cis$PVE_C)*cojo_cis$k)

cojo_cis_unconditional <- cojo_cis[which(cojo_cis$Fstats>=10), ]
cojo_cis_conditional<-cojo_cis[which(cojo_cis$Fstats_C>=10), ]

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

collapsed_df_unconditional <- collapsed_df_unconditional %>%
  dplyr::select(DATASET, TISSUE, SNP, collapsed_Chr, collapsed_bp, collapsed_V2, collapsed_locus_START_END_37,
                collapsed_locus_extended_START_END_37, collapsed_b, collapsed_se, collapsed_mlog10p, collapsed_EA,
                collapsed_NEA, collapsed_MAF, collapsed_freq_geno, collapsed_N, collapsed_PVE, collapsed_k,
                collapsed_Fstats, collapsed_Fstats_multipleMR, collapsed_Entrez_Gene_Name, collapsed_Ensembl_Gene_ID,
                collapsed_TSS, study_id, collapsed_UniProt_ID, collapsed_Target_Name, collapsed_Target_Full_Name,
                FILENAME, Gene.type)


names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Chr"] <- "CHR"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_bp"] <- "POS_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_V2"] <- "POS_38"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_locus_extended_START_END_37"] <- "locus_extended_START_END_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_b"] <- "BETA"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_se"] <- "SE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_mlog10p"] <- "MinusLog10PVAL"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_EA"] <- "EFFECT_ALLELE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_NEA"] <- "OTHER_ALLELE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_MAF"] <- "MAF"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_freq_geno"] <- "EAF"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_N"] <- "SAMPLESIZE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_PVE"] <- "PVE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_k"] <- "k"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Fstats"] <- "Fstats"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "collapsed_Fstats_multipleMR"] <- "Fstats_multipleMR"
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
colnames(collapsed_df_conditional)

collapsed_df_conditional <- collapsed_df_conditional %>%
  dplyr::select(DATASET, TISSUE, SNP, collapsed_Chr, collapsed_bp, collapsed_V2, collapsed_locus_START_END_37,
                collapsed_locus_extended_START_END_37, collapsed_bC, collapsed_bC_se, collapsed_mlog10pC, collapsed_EA, collapsed_NEA, collapsed_MAF,
                collapsed_freq_geno, collapsed_N, collapsed_PVE_C, collapsed_k, collapsed_Fstats_C,
                collapsed_Fstats_C_multipleMR,
                collapsed_Entrez_Gene_Name, collapsed_Ensembl_Gene_ID,
                collapsed_TSS, study_id, collapsed_UniProt_ID,
                collapsed_Target_Name, collapsed_Target_Full_Name, FILENAME, Gene.type)


names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Chr"] <- "CHR"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_bp"] <- "POS_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_V2"] <- "POS_38"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_locus_extended_START_END_37"] <- "locus_extended_START_END_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_bC"] <- "BETA"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_bC_se"] <- "SE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_mlog10pC"] <- "MinusLog10PVAL"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_EA"] <- "EFFECT_ALLELE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_NEA"] <- "OTHER_ALLELE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_MAF"] <- "MAF"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_freq_geno"] <- "EAF"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_N"] <- "SAMPLESIZE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_PVE_C"] <- "PVE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_k"] <- "k"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Fstats_C"] <- "Fstats"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "collapsed_Fstats_C_multipleMR"] <- "Fstats_multipleMR"
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