This folder contains the codes required to download, organise and pre-process ENCODE files. 

Order  of running scripts:  


### Downloading files (FileDownload.sh)
Downloads and organizes .bed files for TF ChIP-seq and ATAC-seq for 15TFs across 5 cell lines. 

Input: metadata.tsv (metadata for of all bed narrowPeak files available for the corresponding experiments on ENCODE)
Output: Each file is saved in the following format: CellLine/Assay/Genome/ENCODEAccessionID_TF.bed.gz. 
Source code: filedownload.py 

### Reference genome annotation file (TSS_referenceFile.py)

### Getting consensus accessibility for each cell line (accessibility_consensus.sh)
Identifies accessible regions in the genome present in all replicates for a given cell line within 10kb of each gene's TSS.

Input: CellLine/ATAC/Genome/ENCODEAccessionID_TF.bed.gz files
Output: CellLine/ATAC/Genome/consensuspeaks_all_targets_10kb.bed files 
Source code: ATAC_consensus.py (Identifies regions present in all replicates)

### Identification of ChIP-seq peaks present in accessible regions (accessible_ChIP_peaks.sh)
Gets set of ChIP-seq peaks that are present in accessible regions for a cell line. Each peak is also assigned to a gene if it is present within 10kbp of the gene's TSS. 

Input: CellLine/ChIP-seq/Genome/ENCODEAccessionID_TF.bed.gz files and CellLine/ATAC/Genome/consensuspeaks_all_targets_10kb.bed
Output: CellLine/ChIP/Genome/TF_consensus_accessible_peaks.bed files and TF_consensus_accessible_peaks_named.bed (Each peak is assigned to genes in vicinity) 


### Fitting the two-state model to the data (TFbindingModelFit.py)
Fits the two-state model to the cumulative positional distribution of TF binding across all genes. Additionally, calculates thresholds for identifying 95% of biologically relevant binding.

Input: CellLine/ATAC/Genome/TF_consensus_accessible_peaks.bed files and TF_consensus_accessible_peaks_named.bed (Each peak is assigned to genes in vicinity) 
Output: processed/Compiled_accessible_target_ChIP_peaks.csv (Cumulative positional distribution of TF binding events)
Processed/ExpDecayConst_params.csv (Parameters for each dataset for two-state model)
Source code: utils.py
