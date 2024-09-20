README MR-instruments CHRIS-INTERVAL METANALYSIS.

20240917 Giulia Pontali, Solène Cadiou, Claudia Giambartolomei @HUMANTECHNOPOLE

GWAS analyses were performed independently using Regenie in our 2 cohorts (INTERVAL study, N = 9251; CHRIS study, N =4194; TOT N = 13445).
Standard-error weighted meta-analysis was performed using METAL software.
Only variants common to the 2 studies were kept in the meta-analysis summary statistics to be analyzed.

We used Locus Breaker function to identify regions of associations (https://www.nature.com/articles/s41467-023-38389-6). 

All analyses including the instrument selection are based on genome build 37. As a last step, to convert positions of the selected instruments to genome build 38, we used the bcftools +liftover tool (https://academic.oup.com/bioinformatics/article/40/2/btae038/7585532).

NOTE 1: (For now missing from the data) Once the MR instruments will be selected for the COLOC step, we will provide the position of the cis region Lifted to build 38.

NOTE 2: We may have multiple cis-instruments corresponding to the same seqID but mapped to a different locus: please treat them as separate signals (not collapsed into one MR model).

To identify cis-instruments for MR, we provide three different files:


1. MR_instruments_best_snps_from_LB_MVP.txt
We first selected as cis instrument the best cis SNP (i.e., SNP with the smallest p-value within a cis locus breaker region) that passes the Bonferroni significancy threshold.
We used the following approach. First, loci from locus breaker results were mapped to cis and trans regions using a cis-trans mapping file. The transcription start site (TSS) information in build 37 from this file was used to define the boundaries of the cis region (+/- 500kb). Any locus overlapping the cis region of the associated protein was defined as cis. We then selected as instrument the SNP having the smallest significant (ie lower than 5.10^-8 divided by our Bonferroni number of effective tests) p-value, if its F-statistics was greater or equal to 10. All top significant SNPs in cis-locus passed the F-statistics filtering, leading to identify 1,799 SNPs in 1,799 cis regions. 

2. MR_instruments_unconditional_analysis_MVP.txt
We used COJO-GCTA to select the conditionally independent SNPs within each locus (extended by +/-100kb) associated to a protein.
We first mapped the SNPs to cis and trans regions using the same cis-trans mapping file as before: any SNP falling in a region defined by +/-500kb from the TSS of the associated protein (as outlined before) was defined as cis.
For the unconditional analysis, we report the unconditional BETA, SE and MinusLog10PVAL. We used unconditional BETA to compute F-statistics (see the formula below).

We identified 4,607 independent SNPs in cis regions. All with F-statistics greater than or equal to 10.

3. MR_instruments_conditional_analysis_MVP.txt
As before, we used COJO-GCTA to select the conditionally independent SNPs within each locus (extended by +/-100kb) associated to a protein.
We first mapped the SNPs to cis and trans regions using a cis-trans mapping file similarly to what has been described in 2.
For the conditional analysis, we report the conditional BETA, SE and MinusLog10PVAL (conditional on all other instruments for that locus). We used conditional BETA to compute F-statistics (see the formula below).
We identified 4,606 independent SNPs in cis regions with F-statistics greater than or equal to 10.

Explanation of column names in the files:

* DATASET = INTERVAL_CHRIS_META
* TISSUE = WholeBlood
* SNP = chr:pos:allele1:allele2 (alleles arranged in alphabetical order)
* CHR
* POS_37 = Position in build 37
* POS_38 = Position in build 38
* locus_START_END_37 = region from locus breaker output
* locus_extended_START_END_37 = region used in fine-mapping (positions from regions from locus breaker extended by +/- 100,000 bp for all regions)
* BETA
* SE
* MinusLog10PVAL
* EFFECT_ALLELE
* OTHER_ALLELE
* MAF
* EAF
* SAMPLESIZE
* PVE = (2*(BETA^2)*MAF*(1-MAF)) / (2*(BETA^2)*MAF*(1-MAF)+(SE^2)*2*SAMPLESIZE*$MAF*(1-MAF))
* k = number of cis conditionally independent SNPs within each locus
* Fstats = PVE*(SAMPLESIZE-1-k)/(1-PVE)*k with k=1
* Fstats_multipleMR = PVE*(SAMPLESIZE-1-k)/(1-PVE)*k (this info is present only in MR_instruments_unconditional_analysis_MVP.txt/MR_instruments_conditional_analysis_MVP.txt, we did not perform any filtering using it)
* GENE_NAME = Please note that when multiple gene names map to the same cis region, the gene names are collapsed (separated by |)
* GENE_ENSEMBL
* TSS_37 = Please note that when multiple gene names map to the same cis region, the TSS_37 positions are collapsed (separated by |)
* SeqID
* UNIPROT
* PROTEIN_NAME
* PROTEIN_LONG_NAME
* FILENAME = NA
* Gene.type = protein_coding
