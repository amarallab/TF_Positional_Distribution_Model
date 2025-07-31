#!/bin/bash
#SBATCH -A b1042 ## <-- EDIT THIS TO BE YOUR ALLOCATION
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --partition  genomics
#SBATCH --job-name=consensusATAC
#SBATCH --mail-type=END,FAIL

module load bedtools
cd TF_Positional_Distribution_Model/data/ENCODE/A549_Perturbation/ATAC-seq



#Run each bed file in a given cell line to get consensus peaks across replicates 
for i in $(ls -d */); do 
	cd $i
	echo $i

	#Output file 
	output_targets="${i%/}_consensuspeaks.bed"
	output="consensuspeaks.bed"

	echo "ATAC-seq initiated"
	gzip -d *bed.gz

	#Sort files
	for bed in $(ls ENC*.bed); do
		out="${bed%.bed}_sorted.bed"
      		sort -V -k1,1 -k2,2 $bed > $out
	done

	#Run multiintersect on all files and select ones with complete consensus
	bedtools multiinter -i ENC*_sorted.bed -header > $output

	#Cleanup files
	rm *sorted.bed
	gzip ENC*.bed

	#Select files with consensus 1
	python ../../../../../scripts/DataProcessing/ENCODE/src/ATAC_consensus.py
        bedtools intersect -a ../../../TF_TSS_10kb_named.bed -b consensuspeaks_all.bed  > $output_targets

	cd ../
done
