rm(list=ls())
library(data.table)

# Upload the mapping file from Eva
mapping_file <- read.delim("/home/giulia.pontali/INTERVAL/mapped_gene_file_GRCh37.txt")
mapping_file <- subset(mapping_file, !duplicated(mapping_file[, c("ID")]))
mapping_file <- mapping_file[ , c("SeqId", "UniProt_ID", "Entrez_Gene_Name", "chromosome", "TSS")]
mapping_file$SeqId <- gsub("-", ".", mapping_file$SeqId)
mapping_file$cis_start <- mapping_file$TSS - 500000
mapping_file$cis_end <- mapping_file$TSS + 500000
ranges <- Map(seq, mapping_file$cis_start, mapping_file$cis_end)
mapping_file$cis_region <- ranges


# Upload COJO's results from Dariush
files <- list.files(path = "", pattern = "\\.txt", full.names = TRUE)

# Loop over the files and read each one
i=1
while(i<=length(files)){
  # Read the gzipped file
  data <- read.delim(files[i])
  seq_number <- sub(".*seq\\.([0-9]+\\.[0-9]+)\\.txt", "\\1", basename(files[i]))
  data$SeqId <- seq_number
  data$Fs <- (data$BETA*data$BETA)/(data$SE*data$SE)
  pass_F_sum_stats <- data[which(data$Fs>10), ] ##check > data[data$F>10, ]

  # Extract cis region:
  merged_df <- merge(pass_F_sum_stats, mapping_file, by = "SeqId", all=F)
  bounds <- strsplit(as.character(merged_df$cis_region), ":")

  merged_df$region <- mapply(function(chr, pos, chromosome, bounds) {
    if (chr == chromosome) {
      range <- as.numeric(bounds)
      return(pos >= range[1] && pos <= range[2])
    } else {
      return(FALSE)
    }
  }, merged_df$CHR, merged_df$POS, merged_df$chromosome, bounds)

  merged_df$region <- ifelse(merged_df$region, "cis", "trans")

  merged_df <- merged_df[which(merged_df$region=="cis"), ]

  if(nrow(merged_df)>1){
    merged_df <- merged_df[,c("SeqId","CHR","POS","SNPID","EA", "NEA", "EAF", "BETA", "SE", "MLOG10P", "Fs", "region")]
    write.table(merged_df, file=file.path("/scratch/giulia.pontali/unconditional_instruments", basename(files[i])),
              row.names = F, quote =F, sep="\t")
  }
  i = i+1
}
