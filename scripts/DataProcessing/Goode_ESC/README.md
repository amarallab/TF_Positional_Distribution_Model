## Description and order of running scripts

###Download raw files (download_fastq.sh)
Download fastq ChIP seq files from GSE69099 (TF ChIP-seq) and GSE69096 (Input files) for different cell types. 
FASTQ files are aligned(bowtie2) against mm39 (indexed bowtie2 file downloaded), duplicates removed (samtools) and saved. 

Input: SRA accession number list (downloaded from SRA run selector).
Output: bam files for each FASTQ file aligned to mm10 and mm39 reference genomes

### Rename aligned reads based on sample ( RenameFiles.py)
Rename SRA accession number to celltype_TF format.

### Peak calling (peak_calling.sh)
Call peaks from bam files using MACS2. 

Input: Aligned reads (.bam files) 
Output: narrowPeaks, summits, peaks.xls. We only use narrow peak files for further analysis. Different p-value thresholds are used for peak calling of mm10 aligned reads. 

### Reference annotation files (TSSreferenceFile.py)
Create a reference file for cds region and +/- 10kb of all TSS for mm10 and mm39. 
Output:TF_TSS_10kb_named_RefGenome.bed and TF_TSS_cds_named_RefGenome.bed

### Identify TF target genes(bedTools_targets.sh)
Find overlapping peaks for TSS in each bed file using bedtools intersect

