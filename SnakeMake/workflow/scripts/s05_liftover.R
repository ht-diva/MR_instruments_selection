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
liftover <- read.table(opt$input_liftover, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)

mapping <- fread(opt$mapping)
lb_liftover_path <- opt$liftover_lb_output

lb$locus_START_END_37 <- paste0(lb$start,"-",lb$end)
lb <- lb[order(lb$CHR, lb$POS), ]

liftover <- liftover[,c(1:3)]
lb <- cbind(liftover, lb)

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
    across(everything(), ~ paste(unique(.), collapse = "|"), .names = "collapsed_{col}"),
    .groups = "drop")

collapsed_df$DATASET="INTERVAL_CHRIS_META_LB"
collapsed_df$TISSUE="WholeBlood"
collapsed_df$FILENAME=NA
collapsed_df$Gene.type = "protein_coding"


collapsed_df <- collapsed_df %>%
  dplyr::select(DATASET, TISSUE, SNPID, collapsed_CHR, collapsed_POS, collapsed_V2, collapsed_locus_START_END_37,
                collapsed_BETA, collapsed_SE, collapsed_MLOG10P, collapsed_EA, collapsed_NEA,
                MAF, collapsed_EAF, collapsed_N, collapsed_PVE, collapsed_Fstats, collapsed_Entrez_Gene_Name, collapsed_Ensembl_Gene_ID,
                collapsed_TSS, phenotype_id, collapsed_UniProt_ID, collapsed_Target_Name, collapsed_Target_Full_Name,
                FILENAME, Gene.type)

names(collapsed_df)[names(collapsed_df) == "collapsed_CHR"] <- "CHR"
names(collapsed_df)[names(collapsed_df) == "collapsed_POS"] <- "POS_37"
names(collapsed_df)[names(collapsed_df) == "collapsed_V2"] <- "POS_38"
names(collapsed_df)[names(collapsed_df) == "collapsed_locus_START_END_37"] <- "locus_START_END_37"
names(collapsed_df)[names(collapsed_df) == "collapsed_BETA"] <- "BETA"
names(collapsed_df)[names(collapsed_df) == "collapsed_SE"] <- "SE"
names(collapsed_df)[names(collapsed_df) == "collapsed_MLOG10P"] <- "MinusLog10PVAL"
names(collapsed_df)[names(collapsed_df) == "collapsed_EA"] <- "EFFECT_ALLELE"
names(collapsed_df)[names(collapsed_df) == "collapsed_NEA"] <- "OTHER_ALLELE"
names(collapsed_df)[names(collapsed_df) == "collapsed_EAF"] <- "EAF"
names(collapsed_df)[names(collapsed_df) == "collapsed_N"] <- "SAMPLESIZE"
names(collapsed_df)[names(collapsed_df) == "collapsed_PVE"] <- "PVE"
names(collapsed_df)[names(collapsed_df) == "collapsed_Fstats"] <- "Fstats"
names(collapsed_df)[names(collapsed_df) == "collapsed_Entrez_Gene_Name"] <- "GENE_NAME"
names(collapsed_df)[names(collapsed_df) == "collapsed_Ensembl_Gene_ID"] <- "GENE_ENSEMBL"
names(collapsed_df)[names(collapsed_df) == "collapsed_TSS"] <- "TSS_37"
names(collapsed_df)[names(collapsed_df) == "collapsed_phenotype_id"] <- "SeqID"
names(collapsed_df)[names(collapsed_df) == "collapsed_UniProt_ID"] <- "UNIPROT"
names(collapsed_df)[names(collapsed_df) == "collapsed_Target_Name"] <- "PROTEIN_NAME"
names(collapsed_df)[names(collapsed_df) == "collapsed_Target_Full_Name"] <- "PROTEIN_LONG_NAME"
names(collapsed_df)[names(collapsed_df) == "collapsed_Gene.type"] <- "Gene.type"

fwrite(collapsed_df, lb_liftover_path)