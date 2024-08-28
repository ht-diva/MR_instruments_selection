suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LB for MR"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--chain_file", default=NULL, help="Chain file to perform liftover"),
  make_option("--liftover_lb_output", default=NULL, help="Output path from MR IVs LB"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
lb <- fread(opt$input)
chain_file <- import.chain(opt$chain_file)
mapping <- fread(opt$mapping)
lb_liftover_path <- opt$liftover_lb_output

# Read the input file
# Prepare the file for liftover
#lb <- read.delim("/home/giulia.pontali/liftover/MR_instruments_best_snps_from_LB.txt", sep = ",")
lb$CHR <- paste0("chr", lb$CHR)

gr <- GRanges(seqnames = lb$CHR,
              ranges = IRanges(start = lb$POS, end = lb$POS),
              strand = "*", names = lb$SNPID)

lifted_data <- liftOver(gr, chain_file)

# Convert to a GRanges object
lifted_granges <- unlist(lifted_data)
lifted_df <- as.data.frame(lifted_granges)

lb <- lb[which(lb$SNPID %in% lifted_df$names), ]

merge_data <- merge(lb, lifted_df, by.x="SNPID", by.y="names", all=F)
merge_data <- merge_data[!duplicated(merge_data), ]
merge_data$CHR <- gsub("chr", "", merge_data$CHR)

mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)

merged <- merge_data %>%
  left_join(mapping, by = c("phenotype_id" = "target"), relationship = "many-to-many") %>%
  filter(CHR == chromosome, (POS >= cis_start & POS <= cis_end))

merged$DATASET="INTERVAL_CHRIS_META_LB"
merged$TISSUE="WholeBlood"
merged$FILENAME=NA
merged$Gene.type = "protein_coding"
merged$MAF <- pmin(merged$EAF, 1-merged$EAF)

merged <- merged %>%
  dplyr::select(DATASET, TISSUE, SNPID, CHR, POS, start, BETA, SE, MLOG10P, EA, NEA, MAF,
         EAF, Fstats, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, phenotype_id, UniProt_ID,
         Target_Name, Target_Full_Name, FILENAME, Gene.type)

names(merged)[names(merged) == "POS"] <- "POS_37"
names(merged)[names(merged) == "start"] <- "POS_38"
names(merged)[names(merged) == "MLOG10P"] <- "MinusLog10PVAL"
names(merged)[names(merged) == "EA"] <- "EFFECT_ALLELE"
names(merged)[names(merged) == "NEA"] <- "OTHER_ALLELE"
names(merged)[names(merged) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged)[names(merged) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged)[names(merged) == "TSS"] <- "TSS_37"
names(merged)[names(merged) == "phenotype_id"] <- "SeqID"
names(merged)[names(merged) == "UniProt_ID"] <- "UNIPROT"
names(merged)[names(merged) == "Target_Name"] <- "PROTEIN_NAME"
names(merged)[names(merged) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"

fwrite(merged, lb_liftover_path)