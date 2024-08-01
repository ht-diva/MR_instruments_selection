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

input_file = snakemake.input.input_file
output_file = snakemake.output.liftover_output_file
unmapped_output_file = snakemake.output.unmapped_output_file

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

# Save the list of SNP IDs before liftover
mysumstats_ids = mysumstats.data[params.snpid]
print(f"Number of SNPs before liftover: {len(mysumstats_ids)}")

# Perform liftover, remove=True remove unmapped variants
mysumstats.liftover(n_cores=3, from_build=params.from_build, to_build=params.to_build, remove=True)

# Save the list of SNP IDs after liftover
mysumstats_ids_liftover = mysumstats.data[params.snpid]
print(f"Number of SNPs after liftover: {len(mysumstats_ids_liftover)}")

unmapped_snps = mysumstats_ids[~mysumstats_ids.isin(mysumstats_ids_liftover)]

unmapped_snps.to_csv(snakemake.output.unmapped_output_file, sep=snakemake.params.sep, index=False)

mysumstats.data.to_csv(snakemake.output.liftover_output_file, sep=snakemake.params.sep, index=False)

