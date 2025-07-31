#!/bin/bash
#SBATCH -A b1042 ## <-- EDIT THIS TO BE YOUR ALLOCATION
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --partition  genomics
#SBATCH --job-name=atac_consensus
#SBATCH --mail-type=END,FAIL

module load bedtools #peak calling

cd TF_Positional_Distribution_Model/data/ENCODE/

for folder in 'A549' 'HepG2' 'GM12878' 'MCF-7' 'K562'; do
	cd $folder
	echo $folder

	### Get ATAC-seq peaks for a cell line
	if [ -d ATAC-seq/ ]; then
		echo "ATAC-seq initiated"
		cd ATAC-seq/
		rm *.bed	#remove all existing processed files
		gzip -d *bed.gz
		for bed in $(ls ENC*.bed); do
			out="${bed%.bed}_sorted.bed"
	      		sort -V -k1,1 -k2,2 $bed > $out
		done
		output="consensuspeaks.bed"
		bedtools multiinter -i ENC*_sorted.bed -header > $output

		rm *sorted.bed
		gzip ENC*.bed
		python ../../../../scripts/DataProcessing/ENCODE/src/ATAC_consensus.py
		output="consensuspeaks_all_targets_10kb.bed"
                bedtools intersect -a ../../TF_TSS_10kb_named.bed -b consensuspeaks_all.bed > $output
		cd ../../
	fi

done
