# pacbioanalysis
Analysis of long read (PacBio) amplicon sequencing data

1. analyse_sam.awk

This is an awk script that accepts a SAM file as an input, and calculates allele/genotype frequencies of different nucleotide and protein variants.
The sequence of the reference DNA is included in the script. 
The script uses a table of genetic code (included in the repository), to translate nucleotide sequences to protein sequences.

2. av_pair_dist.awk

This awk script calculates average pairwise (hamming) distance between protein variants in a population. 
It uses a "haplotype" file (generated by analyse_sam.awk) as an input, that contains a shortened version of the protein sequences in a population. 
A haplotype is stored as a concatenated string of mutated amino acid positions (in ascending order) followed by corresponding mutation. 
For example, 1I157R166R214M is a protein variant that has 4 mutations: M1I, G157R, K166R, K214M. 
An unmutated GFP (ancestral) sequence is represented by an empty string. 


3. GLM_binom.r

This is a R script that identifies mutations that have significantly different frequencies between populations that have evolved under two different conditions (chaperone overexpression and no overexpression).
It takes a table of allele frequencies of single amino acid mutations (rows denote the mutations and columns denote the populations) as input and calculates statistical significance using Generalized Linear Models (binomial family, logit link, likelihood ratio test). 