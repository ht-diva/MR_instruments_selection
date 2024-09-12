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


merged_within_cis <- lb %>%
  left_join(mapping, by = c("phenotype_id" = "SeqId"), relationship = "many-to-many") %>%
  filter(CHR == chromosome) %>%
  mutate(in_cis_range = check_boundaries(start = POS, end = POS, cis_start = cis_start, cis_end = cis_end)) %>%
  filter(in_cis_range == TRUE)

collapsed_df <- merged_within_cis %>%
  # Step 1: Filter rows where in_cis_range == TRUE
  filter(in_cis_range == TRUE) %>%
  # Step 2: Group by phenotype_id and SNPID
  group_by(phenotype_id, SNPID) %>%
  # Step 3: Collapse all columns by concatenating values where necessary
  summarise(
    across(everything(), ~ paste(unique(.), collapse = "|"), .names = "collapsed_{col}")
  ) %>%
  # Step 4: Ungroup to return the result to a normal dataframe
  ungroup()

merged$DATASET="INTERVAL_CHRIS_META_LB"
merged$TISSUE="WholeBlood"
merged$FILENAME=NA
merged$Gene.type = "protein_coding"
merged$MAF <- pmin(merged$EAF, 1-merged$EAF)

merged <- merged %>%
  dplyr::select(DATASET, TISSUE, collapsed_SNPID, collapsed_CHR, collapsed_start, collapsed_end, collapsed_POS, collapsed_V2, collapsed_BETA, collapsed_SE, collapsed_MLOG10P, collapsed_EA, collapsed_NEA,
  MAF, collapsed_EAF, collapsed_N, collapsed_Fstats, collapsed_Entrez_Gene_Name, collapsed_Ensembl_Gene_ID,
  collapsed_TSS, collapsed_phenotype_id, collapsed_UniProt_ID, collapsed_Target_Name, collapsed_Target_Full_Name,
  FILENAME, collapsed_Gene.type)

names(merged)[names(merged) == "collapsed_POS"] <- "POS_37"
names(merged)[names(merged) == "collapsed_V2"] <- "POS_38"
names(merged)[names(merged) == "collapsed_MLOG10P"] <- "MinusLog10PVAL"
names(merged)[names(merged) == "collapsed_EA"] <- "EFFECT_ALLELE"
names(merged)[names(merged) == "collapsed_NEA"] <- "OTHER_ALLELE"
names(merged)[names(merged) == "collapsed_N"] <- "SAMPLESIZE"
names(merged)[names(merged) == "collapsed_Entrez_Gene_Name"] <- "GENE_NAME"
names(merged)[names(merged) == "collapsed_Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(merged)[names(merged) == "collapsed_TSS"] <- "TSS_37"
names(merged)[names(merged) == "collapsed_phenotype_id"] <- "SeqID"
names(merged)[names(merged) == "collapsed_UniProt_ID"] <- "UNIPROT"
names(merged)[names(merged) == "collapsed_Target_Name"] <- "PROTEIN_NAME"
names(merged)[names(merged) == "collapsed_Target_Full_Name"] <- "PROTEIN_LONG_NAME"

fwrite(merged, lb_liftover_path)