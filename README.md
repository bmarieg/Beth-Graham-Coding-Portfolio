# Beth Graham
## Bioinformatics Coding Portfolio

Welcome!

I am using this repository as a coding portfolio (it is not intended as a software package). Each file was included to demonstrate certain skills. I am currently in the process of building it.

Table of Contents:

- "Search and Fetch from Entrez with Biopython":

Jupyter notebook using Biopython to access Entrez databases

- "Rename Fasta Sequences for Easy-to-Read Phylogenetic Trees":

I created this tool to address a personal pet peeve: when you generate a phylogenetic tree from a multiple sequence alignment, the branches are often labeled with accession numbers (as opposed to gene and/or species names), which makes the trees difficult to interpret at a glance. Other times, branches are labeled with obscure scientific species names, which again makes trees difficult to interpret quickly and easily.

This function uses the UniProtKB ID/accession number for each sequence to retrieve the gene name, scientific species name, and the common species name (if available on Entrez). The function then replaces the "description" of each Fasta sequence with an easy-to-read label consisting of the gene and species names. When the newly labeled list of Fasta sequences is used to generate a phylogenetic tree, the tree is much easier to interpret.




Thanks for visiting!

Beth
