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
                         p=params.p,
                         sep=params.sep,
                         verbose=True)

# Perform basic checks
mysumstats.basic_check(verbose=True)

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