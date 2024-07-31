#import gwaslab as gl

# Load the summary statistics - verbose if True print logs
#mysumstats = gl.Sumstats("results/MR_instruments_best_snps_from_LB.txt",
#                         snpid="SNPID",
#                         chrom="CHR",
#                         pos="POS",
#                        ea="EA",
#                        nea="NEA",
#                         neaf="EAF",
#                         beta="BETA",
#                         se="SE",
#                         mlog10p="MLOG10P",
#                         sep=",",
#                         verbose=True)

# Perform basic checks
#mysumstats.basic_check(verbose=False)

# Perform liftover, remove=True remove unmapped variants
#mysumstats.liftover(n_cores=3, from_build="19", to_build="38", remove=True)

# Save the output
#mysumstats.write("MR_instruments_best_snps_from_LB_liftover.txt")
#mysumstats.to_format("results/MR_instruments_best_snps_from_LB_liftover.txt", fmt="gwaslab")

import gwaslab as gl

# Snakemake automatically provides these variables
input_file = snakemake.input[0]
output_file = snakemake.output[0]
unmapped_output_file = snakemake.output[1]  # Add a second output file for unmapped SNPs

params = snakemake.params

# Load the summary statistics
mysumstats = gl.Sumstats(input_file,
                         snpid=params.snpid,
                         chrom=params.chrom,
                         pos=params.pos,
                         ea=params.ea,
                         nea=params.nea,
                         neaf=params.neaf,
                         beta=params.beta,
                         se=params.se,
                         mlog10p=params.mlog10p,
                         sep=params.sep,
                         verbose=True)

# Perform basic checks
mysumstats.basic_check(verbose=False)

# Perform liftover
mysumstats.liftover(n_cores=3, from_build=params.from_build, to_build=params.to_build, remove=True)

# Identify SNPs that do not have corresponding matches (unmapped SNPs)
unmapped_snps = mysumstats.data[mysumstats.data[params.chrom].isnull() | mysumstats.data[params.pos].isnull()]

# Save the unmapped SNPs to a separate file
unmapped_snps.to_csv(unmapped_output_file, sep=params.sep, index=False)

# Save the output for successfully mapped SNPs
mapped_snps = mysumstats.data.dropna(subset=[params.chrom, params.pos])  # Drop SNPs with null chrom or pos
mysumstats.data = mapped_snps
mysumstats.to_format(output_file, fmt="gwaslab")

