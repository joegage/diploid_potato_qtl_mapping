# diploid_potato_qtl_mapping

These are the QTL mapping scripts, genotypic data, and phenotypic data that accompany Marand et al. _Residual heterozygosity and epistatic interactions underlie the complex genetic architecture of yield in diploid potato_. 20XX. doi:xxxxxx

## QTL mapping scripts

QTL mapping in these scripts is performed by fitting multiple linear regression at each recombination bin independently.  The model contains terms for each parental allele (US-W4 and M6), as well as an interaction term.  In bins where one parent is not segregating, that parent and the interaction term are both excluded from the model.

* **1_run_dual_interaction_mapping.R**: This is the primary script that will run QTL mapping for all traits, write the results in text form, and produce plots of significance profiles.
* **map_dual_interactions.R**: contains functions to construct the QTL mapping model (`make_model()`) and to apply the model across all recombination bins for a single phenotype (`dual_mapping()`).
* **run_perm.R**: contains a single function, `runPerm()`, that shuffles the phenotypes and runs QTL mapping.

## Genotypic and phenotypic data

* **W4M6_haplotype_bins_split.bed**: Haplotype designations for each individual at each recombination bin.  Columns 1-3 contain chromosome, haplotype bin start coordinate, and haplotype bin end coordinate. Subsequent pairs of columns contain binary scores for W4 and M6 alleles, respectively. For example, `chr01	0	873946	1	0` corresponds to a bin on chromosome 1 spanning the physical coordinates from 0 to 873946. The individual has W4 haplotype `1` and M6 haplotype `0`. 
* **yield_phenological_traits_allGenos.csv**: contains BLUPs for yield, tuber weight, and tuber number; intercept and slope terms from longitudinal random regression modelling of plant height; and the binary response for growth habit (labelled as 'Height_Group').
