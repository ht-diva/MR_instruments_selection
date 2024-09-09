suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
source("workflow/scripts/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LB for MR"),
  make_option("--vcf_lb_output", default=NULL, help="Output path from vcf"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
df <- fread(opt$input)
vcf_path <- opt$vcf_lb_output

df <- df[order(df$CHR, df$POS), ]

# Function to write VCF file
write_vcf <- function(df, output_filename) {
  # Create a connection to write to a file
  vcf_file <- file(output_filename, "w")

  # Use tryCatch to ensure the file is closed properly
  tryCatch({
    # Write the VCF header
    writeLines("##fileformat=VCFv4.2", vcf_file)
    writeLines("##source=RScript", vcf_file)
    writeLines("##reference=GRCh38", vcf_file)
    writeLines("##INFO=<ID=EAF,Number=1,Type=Float,Description=Effect Allele Frequency>", vcf_file)
    writeLines("##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect Size Estimate>", vcf_file)
    writeLines("##INFO=<ID=SE,Number=1,Type=Float,Description=Standard Error>", vcf_file)
    writeLines("##INFO=<ID=N,Number=1,Type=Integer,Description=Sample Size>", vcf_file)
    writeLines("##INFO=<ID=MLOG10P,Number=1,Type=Float,Description=Negative Log10 P-value>", vcf_file)
    writeLines("##INFO=<ID=phenotype_id,Number=1,Type=String,Description=Phenotype ID>", vcf_file)
    writeLines("##INFO=<ID=cis_or_trans,Number=1,Type=String,Description=Cis or Trans>", vcf_file)
    writeLines("##INFO=<ID=Fstats,Number=1,Type=Float,Description=F-statistics>", vcf_file)
    writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file)

    # Write each row as a VCF entry
    for (i in 1:nrow(df)) {
      # Extract chromosome, position, reference allele (NEA), and alternate allele (EA)
      chrom <- df$CHR[i]
      pos <- df$POS[i]
      id <- df$SNPID[i]
      ref <- df$EA[i]
      alt <- df$NEA[i]
      qual <- "."
      filter <- "."

      # Create the INFO field (customize as needed)
      info <- paste0("EAF=", df$EAF[i], ";BETA=", df$BETA[i], ";SE=", df$SE[i],
                     ";N=", df$N[i], ";MLOG10P=", df$MLOG10P[i],
                     ";phenotype_id=", df$phenotype_id[i],
                     ";cis_or_trans=", df$cis_or_trans[i],
                     ";Fstats=", df$Fstats[i])

      # Write the line to the VCF file
      line <- paste(chrom, pos, id, ref, alt, qual, filter, info, sep = "\t")
      writeLines(line, vcf_file)
    }
  }, finally = {
    # Ensure the file connection is closed
    close(vcf_file)
  })
}

write_vcf(df, vcf_path)

#bgzip output_standard.vcf
#tabix -p vcf output_standard.vcf.gz
#normalization:
#bcftools norm -f human_g1k_v37.fasta -c s -Oz -o data/normalizedtohuman_g1k_v37.vcf.gz data/output_standard.vcf.gz
#bcftools norm -c x -f human_g1k_v37.fasta test_claudia.vcf.gz -o checked.vcf
#bcftools norm -f human_g1k_v37.fasta -c w test_claudia.vcf.gz -Ob -o normalized.vcf.gz (mine does not work)

#bcftools +liftover --no-version -Ou output_standard.vcf.gz -- -s human_g1k_v37.fasta -f hg38.fa -c hg19ToHg38.over.chain.gz > output.lifted.vcf
#bcftools view output.lifted.vcf | less -S
#bcftools view output.lifted.vcf > output.txt

#vcf_data <- read.table("/group//users/giulia.pontali/output.txt", header = F, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)diangelantonio