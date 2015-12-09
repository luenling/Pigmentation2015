# Pigmentation2015

In this repository you can find some of the python and R scripts used for the publication:

**Reconciling differences in Pool-GWAS between populations: a case study of female abdominal pigmentation in Drosophila melanogaster.**

Lukas Endler, Andrea J. Betancourt, Viola Nolte, and Christian Schlötterer; *Genetics* *submitted*

##Python Scripts

### FDR

This directory contains the scripts used for FDR correction. For obtaining an experimental Null distribution, a CMH test is performed on alleles shuffled between control and case according to a beta distribution with a fitting alpha value. The usual workflow is to first obtain the value of alpha giving the P value distribution with the lowest chi squared distance to the observed P value distribution, and afterwards to simulate a Null distribution of 10 times the number of tested variants with this alpha value to estimate false detection rates. 

For this run the script `null_dist3_alt.py` on a CMH or sync file giving the the replicates to compare in the CMH test with with varying alpha values to get the chi-square distance between observed and shuffled P values. For this you can either use alpha == beta for the beta-binomial sampling or give your population coverages to get an estimate for a
beta (only if the coverages are very different). Create qq plots for the best distributions to see how great they look.

     python null_dist3_alt.py --alpha 25 --in infile.sync --chi2 --reps "1,2,3,1,2,3" >> fdr.log

This creates creates a file with P values called sync/cmh file + alpha/beta values + "`_nullP.out`" and adds to a file called after the used sync/cmh file + "`chisquare.out`", in which each line gives the alpha value, the beta values for each replicate, and the chi-square distance between the observed and shuffled P values. For the alpha value giving the lowest chi squared value, rerun the script with the option --10x to create a file with 10 times the P values.

     python null_dist3_alt.py --alpha 25 --in infile.sync --chi2 --10x --reps "1,2,3,1,2,3" >> fdr.log

With the resulting file with the extension `_nullP.out_10x` you can use the script `add_fdr_q_values.py` to calculate
and add the q values for your synchronised file containing the P value in the last column (tab delimited: Chromosome, Baseposition, Reference allele, ... , P-value ). This results in a file with Chromosome, Baseposition, Reference allele,  P-value, and q value, sorted by increasing q values.

    python add_fdr_q_values.py -s infile.pvalues -n infile_nullP.out_10x -o outfile.qvalues

### `get_inversion_freqs_nosyncparser.py` and `Dmel_inversions_tab.SNPmarkers`

Using a tab delimited file with the alleles fixed in an inversion ( see: `Dmel_inversions_tab.SNPmarkers`) together with a sync file, `get_inversion_freqs_nosyncparser.py` allows to estitmate the frequency of an inversion in the pooled samples. It gives a file with the positions of alleles fixed in the inversions, the reference and inverted allele, and the fractions of inverted allele found in the different pools in the sync file. For estimating the inversion frequencies we commonly use the median of fractions for each inversion. 

    python get_inversion_freqs_nosyncparser.py --inv Dmel_inversions_tab.SNPmarkers --infile infile.sync --header > inversion_fractions.sync 

### `cmh_test.py` and `cmh_test_indels.py`

Both scripts create P values, odds ratio and conf intervals for a CMH test on a sync file or an indel count sync file.
For a SNP synvc file `cmh_test.py` always compares the overal most common (major) allele against the overall second most common (minor). The third column in the output file gives the major and minor alleles used for the test.


The indel count sync file should have a format akin to a tab delimited sync file (Chromosome
position allele1 allele2,3,4... (count_a1:count_a2:count:...)*n) and the odds ratio is always given for the comparison a1:a2.



### Sync_parser.py

This module contains the class `SyncLineParser` for reading sync files as produced by [`PoPoolation2`](http://sourceforge.net/projects/popoolation2/). Integrated into all python scripts, so not necessary as a separate module.




When given a 

FDR                                 cmh_test.py                         extract_haplos_from_consensus.py    get_inversion_freqs_nosyncparser.py sample_SNPs_for_snp_overlap.py      sample_r2_for_snp_pairs.py
README.md                           count_ancestral_polymorphism.py     get_ancestral_allele_from_sim.py    get_light_dark_alleles.py           sample_SNPs_rank_lOR_sign.py        simulate_qtl_uniform.R
Sync_parser.py                      create_list_of_LL_DD_haps.py        get_haplotypes.py                   sample_SNPs_for_gene_overlap.py     sample_ancestral_polymorph.py
