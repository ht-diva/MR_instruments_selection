suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of Cojo"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--chain_file", default=NULL, help="Chain file to perform liftover"),
  make_option("--conditional_output", default=NULL, help="Output path and name for list of instruments from Cojo"),
  make_option("--unconditional_output", default=NULL, help="Output path and name for list of instruments from Cojo unconditional"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cojo<-fread(opt$input)
chain_file <- import.chain(opt$chain_file)
conditional_path<-opt$conditional_output
unconditional_path<-opt$unconditional_output
mapping<-fread(opt$mapping)

##mapping file
# mapping<-fread("/home/solene.cadiou/basic_GWAS_protein/meta_results/MR/MR_instruments_selection/mapped_gene_file_GRCh37_21052025.txt")
mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)
##################
###cis mapping of Cojo file
for (i in 1:nrow(cojo)){
  cojo$cis_or_trans[i]<-map(cojo$study_id[i],cojo$Chr[i],cojo$bp[i],cojo$bp[i],mapping)
}
##filtering only cis
cojo_cis<-cojo[which(cojo$cis_or_trans=="cis"),]
##select unconditional beta se
cojo_cis$Fstats_C<-((cojo_cis$bC^2)/(cojo_cis$bC_se^2))
cojo_cis<-cojo_cis[which(cojo_cis$Fstats_C>=10), ]
cojo_cis$Fstats<-((cojo_cis$b^2)/(cojo_cis$se^2))

cojo_cis$Chr <- paste0("chr", cojo_cis$Chr)

gr <- GRanges(seqnames = cojo_cis$Chr,
              ranges = IRanges(start = cojo_cis$bp, end = cojo_cis$bp),
              strand = "*", names = cojo_cis$SNP)

lifted_data <- liftOver(gr, chain_file)

lifted_granges <- unlist(lifted_data)
lifted_df <- as.data.frame(lifted_granges)

cojo_cis <- cojo_cis[which(cojo_cis$SNP %in% lifted_df$names), ]

merge_data <- merge(cojo_cis, lifted_df, by.x="SNP", by.y="names", all=F)
merge_data <- merge_data[!duplicated(merge_data), ]
merge_data$Chr <- gsub("chr", "", merge_data$Chr)

split_variants <- strsplit(merge_data$SNP, ":")

# Extract EA (Effect Allele) and NEA (Non-Effect Allele)
merge_data$EA <- sapply(split_variants, function(x) x[3])
merge_data$NEA <- sapply(split_variants, function(x) x[4])

merge_data$MAF <- pmin(merge_data$freq, 1-merge_data$freq)
merge_data$DATASET <- "INTERVAL_CHRIS_META_COJO"
merge_data$TISSUE <- "WholeBlood"
merge_data$FILENAME <- NA
merge_data$Gene.type <- "protein_coding"

merged <- merge_data %>%
  left_join(mapping, by = c("study_id" = "target"), relationship = "many-to-many") %>%
  filter(Chr == chromosome, (bp >= cis_start & bp <= cis_end))
merged <- merged[!duplicated(merged), ]

merged_conditional <- merged %>%
  dplyr::select(DATASET, TISSUE, SNP, Chr, bp, start, bC, bC_se, mlog10pC, EA, NEA, MAF,
                freq, Fstats_C, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, study_id, UniProt_ID,
                Target_Name, Target_Full_Name, FILENAME, Gene.type)

names(merged_conditional)[names(merged_conditional) == "SNP"] <- "SNPID"
names(merged_conditional)[names(merged_conditional) == "Chr"] <- "CHR"
names(merged_conditional)[names(merged_conditional) == "bp"] <- "POS_37"
names(merged_conditional)[names(merged_conditional) == "start"] <- "POS_38"
names(merged_conditional)[names(merged_conditional) == "bC"] <- "BETA_conditional"
names(merged_conditional)[names(merged_conditional) == "bC_se"] <- "SE_conditional"
names(merged_conditional)[names(merged_conditional) == "mlog10pC"] <- "MinusLog10PVAL_conditional"
names(merged_conditional)[names(merged_conditional) == "EA"] <- "EFFECT_ALLELE"
names(merged_conditional)[names(merged_conditional) == "NEA"] <- "OTHER_ALLELE"
names(merged_conditional)[names(merged_conditional) == "freq"] <- "EAF"
names(merged_conditional)[names(merged_conditional) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged_conditional)[names(merged_conditional) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged_conditional)[names(merged_conditional) == "TSS"] <- "TSS_37"
names(merged_conditional)[names(merged_conditional) == "study_id"] <- "SeqID"
names(merged_conditional)[names(merged_conditional) == "UniProt_ID"] <- "UNIPROT"
names(merged_conditional)[names(merged_conditional) == "Target_Name"] <- "PROTEIN_NAME"
names(merged_conditional)[names(merged_conditional) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"


cojo_conditional<-merged_conditional

merged_unconditional <- merged %>%
  dplyr::select(DATASET, TISSUE, SNP, Chr, bp, start, b, se, MLOG10P, EA, NEA, MAF,
                freq, Fstats, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, study_id, UniProt_ID,
                Target_Name, Target_Full_Name, FILENAME, Gene.type)

names(merged_unconditional)[names(merged_unconditional) == "SNP"] <- "SNPID"
names(merged_unconditional)[names(merged_unconditional) == "Chr"] <- "CHR"
names(merged_unconditional)[names(merged_unconditional) == "bp"] <- "POS_37"
names(merged_unconditional)[names(merged_unconditional) == "start"] <- "POS_38"
names(merged_unconditional)[names(merged_unconditional) == "b"] <- "BETA"
names(merged_unconditional)[names(merged_unconditional) == "se"] <- "SE"
names(merged_unconditional)[names(merged_unconditional) == "MLOG10P"] <- "MinusLog10PVAL"
names(merged_unconditional)[names(merged_unconditional) == "EA"] <- "EFFECT_ALLELE"
names(merged_unconditional)[names(merged_unconditional) == "NEA"] <- "OTHER_ALLELE"
names(merged_unconditional)[names(merged_unconditional) == "freq"] <- "EAF"
names(merged_unconditional)[names(merged_unconditional) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged_unconditional)[names(merged_unconditional) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged_unconditional)[names(merged_unconditional) == "TSS"] <- "TSS_37"
names(merged_unconditional)[names(merged_unconditional) == "study_id"] <- "SeqID"
names(merged_unconditional)[names(merged_unconditional) == "UniProt_ID"] <- "UNIPROT"
names(merged_unconditional)[names(merged_unconditional) == "Target_Name"] <- "PROTEIN_NAME"
names(merged_unconditional)[names(merged_unconditional) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

cojo_unconditional<-merged_unconditional

fwrite(cojo_conditional, conditional_path)
fwrite(cojo_unconditional,unconditional_path)
fwrite(cojo, cojo_path)