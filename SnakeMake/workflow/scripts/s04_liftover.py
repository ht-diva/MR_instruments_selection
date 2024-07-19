import gwaslab as gl

# Load the summary statistics - verbose if True print logs
mysumstats = gl.Sumstats("MR_instruments_best_snps_from_LB.txt",
                         snpid="SNP",
                         chrom="CHR",
                         pos="POS",
                         ea="ALT",
                         nea="REF",
                         neaf="Frq",
                         beta="BETA",
                         se="SE",
                         p="P",
                         verbose=True)

# Perform basic checks
mysumstats.basic_check(verbose=False)

# Perform liftover, remove=True remove unmapped variants
mysumstats.liftover(n_cores=3, from_build="19", to_build="38", remove=True)

# Save the output
mysumstats.write("MR_instruments_best_snps_from_LB_liftover.txt")