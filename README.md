# MR_instruments_selection

## Strategy:

1. best SNP from pvalue in locus breaker region overlapping cis region and passing bonferroni threshold
     * input: locusbreaker results + sumstats
     * parameters/other imputs: Bonferroni threshold, cis trans mapping file
     * steps:
         * cis/trans map the locus breaker results
         * filter the cis
         * compute the f stats for the best snps from LB results
         * if F stats >=10, select the best snps
         * if not:\
               1. extract the locusbreaker region from summary statistics\
               2. select only the snp passing Bonferroni threshold\
               3. compute fstats and filter out those with Fstats <=10\
               4. select the most significant snp if any remaining
     * output: list of instruments (<=1 per target)

2. insturments from conditional using conditional beta and SE estimates
     * input: conditional sumstats
     * parameters/other inputs: cis trans mapping file
     * steps:
         * mapping: output is sumstats with additionaval cis trans variable
         * filtering cis\
         (Alternative strategy: take end and start (of cis region defined by TSS) from mapping file and extract the cis region)\
         * calculate fstats from conditional beta and SE
         * filter out fstats <=10
         * save remaining snps (with beta and SE pvalue)

3. insturments from conditional using unconditional beta and SE estimates
     * input: instruments defined on conditional sumstats, unconditional sumstats
     * parameters/other imputs:
     * steps
         * keep only the snps in the previous list of instruments
         * save remaining snps (with beta and SE pvalue)
