import gwaslab as gl

# Load the summary statistics - verbose if True print logs
mysumstats = gl.Sumstats("results/MR_instruments_best_snps_from_LB.txt",
                         snpid="SNPID",
                         chrom="CHR",
                         pos="POS",
                         ea="EAF",
                         nea="NEA",
                         neaf="EAF",
                         beta="BETA",
                         se="SE",
                         p="MLOG10P",
                         sep=",",
                         verbose=True) #what this verbose do?

# See how to add how many SNPs we have and how many snps we loose after harmonization.
# Perform basic checks
mysumstats.basic_check(verbose=False)

# Perform liftover, remove=True remove unmapped variants
mysumstats.liftover(n_cores=3, from_build="19", to_build="38", remove=True)

# Save the output
mysumstats.write("results/MR_instruments_best_snps_from_LB_liftover.txt")