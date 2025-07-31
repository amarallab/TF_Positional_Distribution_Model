#!/bin/bash
#SBATCH -A b1042 ## <-- EDIT THIS TO BE YOUR ALLOCATION
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --partition  genomics
#SBATCH --job-name=chip_consensus
#SBATCH --mail-type=END,FAIL

module load bedtools #peak calling

cd TF_Positional_Distribution_Model/data/ENCODE/

for folder in 'A549' 'HepG2' 'GM12878' 'MCF-7' 'K562'; do
	cd $folder
	echo $folder
	
	#Sort ATAC-seq file for multiinter
	sort -k 1,1 -k2,2n ATAC-seq/consensuspeaks_all_targets_10kb.bed > ATAC-seq/consensuspeaks_all_targets_10kb_sorted.bed
	
	cd ChIP-seq
	
	#For each TF get set of consensus accessible peaks
	for tf in 'CEBPB' 'CTCF' 'ELK1' 'EP300' 'ESRRA' 'HDAC2' 'MAFK' 'MAX' 'MAZ' 'MYC' 'RAD21' 'RCOR1' 'REST' 'RFX5' 'SIN3A'; do
		
		#Unzip each raw bed file
		gzip -d *$tf*

		#Sort file
		for bed in $(ls *$tf*); do
			out="${bed%%.bed}_sorted.bed"
			sort -k 1,1 -k2,2n $bed  > $out
		done

		#Run multiinter
		output="${tf}_consensus_accessible_peaks.bed"
		bedtools multiinter -i *$tf*sorted.bed ../ATAC-seq/consensuspeaks_all_targets_10kb_sorted.bed -header > $output
		
		#Get corresponding gene for each peak
		output_named="${output%%.bed}_named.bed"
		bedtools intersect -a ../../TF_TSS_10kb_named.bed -b $output > $output_named
		
		#Remove intermediate sort files
		rm *sorted.bed
		gzip ENC*$tf*
	done
	cd ../../
done
