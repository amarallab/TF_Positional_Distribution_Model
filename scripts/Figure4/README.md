This folder contains the scripts needed to generate Figure 4 of the paper.

## Downloading the data
To download the data from ENCODE, run the following line of code in the 
TF_Positional_Distribution_Model/data/ENCODE/RNA-seq/ directory

xargs -n 1 curl -O -L < files.txt

Output: ENCODEAccessionID.tsv (each accession ID is 1 RNA-seq sample)

## Compiling RNA-seq replicates (RNAseq_TPM_compile.py)
This script combines the sequencing data for all replicates, converts ENSEMBL IDs to gene symbols and calculates TPM counts
for each gene. 
Input: ENCODEAccessionID.tsv
Output: Compiled_RNAseq.csv

## Generating Figure 4 (Figure4.ipynb)
This notebook generates Figure 4 in the paper. For each TF, biologically relevant and spurious target genes are classified.
Next, for each pair of TF-target gene, the correlation in mRNA expression is estimated (Spearman's correlation coefficient). 



