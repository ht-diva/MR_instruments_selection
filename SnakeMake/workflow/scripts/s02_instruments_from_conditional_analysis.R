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
cojo$locus_split_start <- as.numeric(cojo$locus_split_start)
cojo$locus_split_end <- sapply(locus_split, function(x) x[2])
cojo$locus_split_end <- as.numeric(cojo$locus_split_end)
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
  cojo$cis_or_trans[i]<-map(cojo$study_id[i],cojo$Chr[i],cojo$locus_split_start[i],cojo$locus_split_end[i],mapping)
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

collapsed_df_unconditional <- cojo_cis_unconditional %>%
  left_join(mapping %>% distinct(target, chromosome, .keep_all = TRUE), by = c("study_id" = "target"), relationship = "many-to-many") %>%
  filter(Chr == chromosome,
  locus_split_start <= cis_end & locus_split_end >= cis_start
  )

collapsed_df_unconditional$DATASET="INTERVAL_CHRIS_META_COJO_unconditional"
collapsed_df_unconditional$TISSUE="WholeBlood"
collapsed_df_unconditional$FILENAME=NA
collapsed_df_unconditional$Gene.type = "protein_coding"

collapsed_df_unconditional <- collapsed_df_unconditional %>%
  dplyr::select(DATASET, TISSUE, SNP, Chr, bp, V2, locus_START_END_37,
                locus_extended_START_END_37, b, se, mlog10p, EA,
                NEA, MAF, freq_geno, N, PVE, k,
                Fstats, Fstats_multipleMR, Entrez_Gene_ID,
                TSS, study_id, UniProt_ID, Target_Name, Target_Full_Name,
                FILENAME, Gene.type)


names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "Chr"] <- "CHR"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "bp"] <- "POS_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "V2"] <- "POS_38"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "locus_extended_START_END_37"] <- "locus_extended_START_END_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "b"] <- "BETA"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "se"] <- "SE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "mlog10p"] <- "MinusLog10PVAL"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "EA"] <- "EFFECT_ALLELE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "NEA"] <- "OTHER_ALLELE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "MAF"] <- "MAF"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "freq_geno"] <- "EAF"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "N"] <- "SAMPLESIZE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "PVE"] <- "PVE"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "k"] <- "k"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "Fstats"] <- "Fstats"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "Fstats_multipleMR"] <- "Fstats_multipleMR"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "Entrez_Gene_ID"] <- "GENE_NAME"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "TSS"] <- "TSS_37"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "study_id"] <- "SeqID"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "UniProt_ID"] <- "UNIPROT"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "Target_Name"] <- "PROTEIN_NAME"
names(collapsed_df_unconditional)[names(collapsed_df_unconditional) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

cojo_unconditional<-collapsed_df_unconditional

collapsed_df_conditional <- cojo_cis_conditional %>%
  left_join(mapping %>% distinct(target, chromosome, .keep_all = TRUE), by = c("study_id" = "target"), relationship = "many-to-many") %>%
  filter(Chr == chromosome,
  locus_split_start <= cis_end & locus_split_end >= cis_start
  )

collapsed_df_conditional$DATASET="INTERVAL_CHRIS_META_COJO_conditional"
collapsed_df_conditional$TISSUE="WholeBlood"
collapsed_df_conditional$FILENAME=NA
collapsed_df_conditional$Gene.type = "protein_coding"
colnames(collapsed_df_conditional)

collapsed_df_conditional <- collapsed_df_conditional %>%
  dplyr::select(DATASET, TISSUE, SNP, Chr, bp, V2, locus_START_END_37,
                locus_extended_START_END_37, bC, bC_se, mlog10pC, EA, NEA, MAF,
                freq_geno, N, PVE_C, k, Fstats_C,
                Fstats_C_multipleMR,
                Entrez_Gene_ID,
                TSS, study_id, UniProt_ID,
                Target_Name, Target_Full_Name, FILENAME, Gene.type)


names(collapsed_df_conditional)[names(collapsed_df_conditional) == "Chr"] <- "CHR"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "bp"] <- "POS_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "V2"] <- "POS_38"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "locus_extended_START_END_37"] <- "locus_extended_START_END_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "bC"] <- "BETA"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "bC_se"] <- "SE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "mlog10pC"] <- "MinusLog10PVAL"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "EA"] <- "EFFECT_ALLELE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "NEA"] <- "OTHER_ALLELE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "MAF"] <- "MAF"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "freq_geno"] <- "EAF"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "N"] <- "SAMPLESIZE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "PVE_C"] <- "PVE"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "k"] <- "k"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "Fstats_C"] <- "Fstats"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "Fstats_C_multipleMR"] <- "Fstats_multipleMR"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "Entrez_Gene_ID"] <- "GENE_NAME"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "TSS"] <- "TSS_37"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "study_id"] <- "SeqID"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "UniProt_ID"] <- "UNIPROT"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "Target_Name"] <- "PROTEIN_NAME"
names(collapsed_df_conditional)[names(collapsed_df_conditional) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

cojo_conditional<-collapsed_df_conditional

fwrite(cojo_conditional, conditional_path)
fwrite(cojo_unconditional,unconditional_path)