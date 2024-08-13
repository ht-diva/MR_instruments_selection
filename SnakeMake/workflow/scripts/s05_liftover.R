#install.packages("rtracklayer")
rm(list=ls())
library(rtracklayer)

# Read the input file
# Prepare the file for liftover
a <- read.delim("/home/giulia.pontali/liftover/MR_instruments_best_snps_from_LB.txt", sep = ",")
b <- a[, c(1, 2, 2, 3)]
b$POS.1 <- b$POS.1 + 1
b$CHR <- paste0("chr", b$CHR)

#write.table(b, "/home/giulia.pontali/liftover_file.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

chain_file <- import.chain("/home/giulia.pontali/liftover/hg19ToHg38.over.chain")
# Read BED file
bed_data <- import("liftover_file.bed", format = "BED")
lifted_data <- liftOver(bed_data, chain_file)

# Convert to a GRanges object
lifted_granges <- unlist(lifted_data)
lifted_df <- as.data.frame(lifted_granges)

condition <- a$SNPID == lifted_df$name
result_list = list()

for (i in 1:nrow(a)) {
  if (condition[i]) {
    # Full merge when condition is TRUE
      merged_row = cbind(a[i, ], lifted_df[i, !(names(lifted_df) %in% names(a))])
      } else {
        # Partial merge on common columns when condition is FALSE
          common_cols = intersect(names(a), names(lifted_df))
          merged_row = merge(a[i, common_cols, drop = FALSE], lifted_df[i, common_cols, drop = FALSE])
          }

    # Add the merged row to the list
    result_list[[i]] = merged_row
}

merged_df = do.call(rbind, result_list)

mapping <- read.delim("/home/giulia.pontali/INTERVAL/mapped_gene_file_GRCh37.txt")
mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)

merged <- merged_df %>%
  left_join(mapping, by = c("phenotype_id" = "target"), relationship = "many-to-many")

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
names(merged)[names(merged) == "NEA"] <- "OTHER_ALLELE"
names(merged)[names(merged) == "Entrez_Gene_Name"] <- "GENE_NAME"
names(merged)[names(merged) == "Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged)[names(merged) == "TSS"] <- "TSS_37"
names(merged)[names(merged) == "phenotype_id"] <- "SeqID"
names(merged)[names(merged) == "UniProt_ID"] <- "UNIPROT"
names(merged)[names(merged) == "Target_Name"] <- "PROTEIN_NAME"
names(merged)[names(merged) == "Target_Full_Name"] <- "PROTEIN_LONG_NAME"
