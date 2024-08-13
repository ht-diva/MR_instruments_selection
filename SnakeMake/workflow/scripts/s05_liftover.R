suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LB for MR"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--chain_file", default=NULL, help="Chain file to perform liftover")
  make_option("--liftover_lb_output", default=NULL, help="Output path from MR IVs LB"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
lb <- fread(opt$input)
chain_file <-fread(opt$chain_file)
mapping <- fread(opt$mapping)
lb_liftover_path <- opt$liftover_lb_output

# Read the input file
# Prepare the file for liftover
#lb <- read.delim("/home/giulia.pontali/liftover/MR_instruments_best_snps_from_LB.txt", sep = ",")
for_liftover <- lb[, c(1, 2, 2, 3)]
for_liftover$POS.1 <- for_liftover$POS.1 + 1
for_liftover$CHR <- paste0("chr", for_liftover$CHR)

gr <- GRanges(seqnames = for_liftover$CHR,
              ranges = IRanges(start = for_liftover$POS.1, end = for_liftover$POS.2),
              strand = "*")

lifted_data <- liftOver(gr, chain_file)

# Convert to a GRanges object
lifted_granges <- unlist(lifted_data)
lifted_df <- as.data.frame(lifted_granges)

condition <- lb$SNPID == lifted_df$name
result_list = list()

for (i in 1:nrow(lb)) {
  if (condition[i]) {
    # Full merge when condition is TRUE
      merged_row = cbind(lb[i, ], lifted_df[i, !(names(lifted_df) %in% names(lb))])
      } else {
        # Partial merge on common columns when condition is FALSE
          common_cols = intersect(names(lb), names(lifted_df))
          merged_row = merge(lb[i, common_cols, drop = FALSE], lifted_df[i, common_cols, drop = FALSE])
          }

    # Add the merged row to the list
    result_list[[i]] = merged_row
}

merged_df = do.call(rbind, result_list)

mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")

merged <- merged_df %>%
  left_join(mapping, by = c("phenotype_id" = "target"), relationship = "many-to-many") %>%
  filter(CHR == chromosome)

merged$DATASET="INTERVAL_CHRIS_META_LB"
merged$TISSUE="WholeBlood"
merged$FILENAME=NA
merged$Gene.type = "protein_coding"

merged <- merged %>%
  select(DATASET, TISSUE, SNPID, CHR, POS, start, BETA, SE, MLOG10P, EA, NEA,
         EAF, Entrez_Gene_Name, Ensembl_Gene_ID, TSS, phenotype_id, UniProt_ID,
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