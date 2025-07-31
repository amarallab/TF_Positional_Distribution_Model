This folder contains the scripts needed to generate Figure 5 of the paper.

## Downloading the data (FileDownload.py)
Files associated with 100nm dexamethasone (12h) treatment in A549 cells are downloaded. This includes ChIP-seq data for 7TFs, ATAC-seq and RNA-seq 
data for control and dexamethasone treated cells.
Output: Assay/Treatment/TF/ENCODEAccessionID.bed (ChIP-seq file) or Assay/Treatment/ENCODEAccessionID.bed (ATAC-seq) or.tsv (RNA-seq)

### Getting consensus accessibility for each cell line (accessibility_consensus.sh)
Identifies accessible regions in the genome present in all replicates for a given cell line within 10kb of each gene's TSS for each condition.

Input: A549_Perturbation/ATAC/Treatment/ENCODEAccessionID_TF.bed.gz files
Output: A549_Perturbation/ATAC/Treatment/consensuspeaks_all_targets_10kb.bed files 
Source code: ATAC_consensus.py (Identifies regions present in all replicates)

### Identification of ChIP-seq peaks present in accessible regions (accessible_ChIP_peaks.sh)
Gets set of ChIP-seq peaks that are present in accessible regions for a cell line. Each peak is also assigned to a gene if it is present within 10kbp of the gene's 
TSS. 

Input: A549_Perturbation/ChIP-seq/Treatment/TF/ENCODEAccessionID_TF.bed.gz files and A549_Perturbation/ATAC/Treatment/consensuspeaks_all_targets_10kb.bed
Output: A549_Perturbation/ChIP-seq/Treatment/TF/TF_consensus_accessible_peaks.bed files and TF_consensus_accessible_peaks_named.bed (Each peak is assigned to genes in 
vicinity) 


### Fitting the two-state model to the data (TFbindingModelFit.py)
Fits the two-state model to the cumulative positional distribution of TF binding across all genes. Additionally, calculates thresholds for identifying 95% of 
biologically relevant binding.

Input: A549_Perturbation/ChIP-seq/Treatment/TF/TF_consensus_accessible_peaks.bed
Output: processed/Compiled_accessible_target_ChIP_peaks.csv (Cumulative positional distribution of TF binding events)
Processed/ExpDecayConst_params.csv (Parameters for each dataset for two-state model)
Source code: utils.py


## RNA-seq data processing (RNAseq_compile_TPM.py)
This script combines the sequencing data for all replicates, converts ENSEMBL IDs to gene symbols and calculates TPM counts
for each gene. 
Input: ENCODEAccessionID.tsv
Output: A549_dex12h_RNAseq_counts.csv,A549_dex12h_RNAseq_TPM.csv


## Generating Figure 5 (Figure5.ipynb)
This notebook generates Figure 5 in the paper. 



