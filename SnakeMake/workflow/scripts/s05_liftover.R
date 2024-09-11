suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LB for MR"),
  make_option("--input_liftover", default=NULL, help="Path and file name of liftover of LB"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--liftover_lb_output", default=NULL, help="Output path from MR IVs LB"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
lb <- fread(opt$input)
liftover <- fread(opt$input_liftover)
mapping <- fread(opt$mapping)
lb_liftover_path <- opt$liftover_lb_output

lb <- lb[order(lb$CHR, lb$POS), ]
liftover <- liftover[,c(1:3)]

data <- cbind(liftover, lb)

mapping <- subset(mapping, !duplicated(mapping[, c("ID")]))
mapping$SeqId <- gsub("-", ".", mapping$SeqId)
mapping$SeqId <- paste0("seq.", mapping$SeqId )
mapping$cis_start <- mapping$TSS - 500000
mapping$cis_end <- mapping$TSS + 500000


merged <- data %>%
  left_join(mapping, by = c("phenotype_id" = "SeqId"), relationship = "many-to-many") %>%
  filter(CHR == chromosome, (POS >= cis_start & POS <= cis_end))

merged$DATASET="INTERVAL_CHRIS_META_LB"
merged$TISSUE="WholeBlood"
merged$FILENAME=NA
merged$Gene.type = "protein_coding"
merged$MAF <- pmin(merged$EAF, 1-merged$EAF)

merged <- merged %>%
  dplyr::select(DATASET, TISSUE, SNPID, CHR, POS, V2, BETA, SE, MLOG10P, EA, NEA, MAF,
                EAF, N, Fstats, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, phenotype_id, UniProt_ID,
                Target_Name, Target_Full_Name, FILENAME, Gene.type)

names(merged)[names(merged) == "POS"] <- "POS_37"
names(merged)[names(merged) == "V2"] <- "POS_38"
names(merged)[names(merged) == "MLOG10P"] <- "MinusLog10PVAL"
names(merged)[names(merged) == "EA"] <- "EFFECT_ALLELE"
names(merged)[names(merged) == "NEA"] <- "OTHER_ALLELE"
names(merged)[names(merged) == "N"] <- "SAMPLESIZE"
names(merged)[names(merged) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged)[names(merged) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged)[names(merged) == "TSS"] <- "TSS_37"
names(merged)[names(merged) == "phenotype_id"] <- "SeqID"
names(merged)[names(merged) == "UniProt_ID"] <- "UNIPROT"
names(merged)[names(merged) == "Target_Name"] <- "PROTEIN_NAME"
names(merged)[names(merged) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

fwrite(merged, lb_liftover_path)