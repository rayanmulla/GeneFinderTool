# GeneFinderTool
A gene finder tool made in python

gene_finder.py outputs open reading frames (ORFs) by scanning a FASTA file (.fna)

rosalind_72.py outputs utilizes the code in gene_finder.py to solve Rosalind problem number 72 outputing protein strings that can be translated from the ORFs

Use the coding files (find_orfs.sh) and (process_orfs.py) inside the folder data to output ORFs in all 14 files of genomes inside the folder with two filters applied: minimum length (default minimum ORF length = 100 codons) and RBS detection that searches for Shine-Dalgarno sequence upstrem the starting codon. 

