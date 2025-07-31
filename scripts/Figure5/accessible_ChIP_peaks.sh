#!/bin/bash
#SBATCH -A b1042 ## <-- EDIT THIS TO BE YOUR ALLOCATION
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --partition  genomics
#SBATCH --job-name=consensus
#SBATCH --mail-type=END,FAIL

module load bedtools #peak calling
cd TF_Positional_Distribution_Model/data/ENCODE/A549_Perturbation/ChIP-seq

#Sort ATAC seq files
out="../ATAC-seq/control/control_consensuspeaks_sorted.bed"
sort -V -k1,1 -k2,2 "../ATAC-seq/control/control_consensuspeaks.bed" > $out

out="../ATAC-seq/dexamethasone/dexamethasone_consensuspeaks_sorted.bed"
sort -V -k1,1 -k2,2 "../ATAC-seq/dexamethasone/dexamethasone_consensuspeaks.bed" > $out

#Run each bed file in a given TF,treatment to get consensus peaks across replicates
for i in $(ls -d */); do 
	cd $i

	for tf in $(ls -d */); do
		cd $tf

		rm *.bed
		#Output files for identifying accessible TF binding sites
		output_dex="${i%/}_consensuspeaks_accessible_dexamethasone.bed"
                output_control="${i%/}_consensuspeaks_accessible_control.bed"

		#Unzip and sort files
		gzip -d *.gz
		for bed_file in *.bed; do
        		out="${bed_file%.bed}_sorted.bed"
        	        sort -V -k1,1 -k2,2 $bed_file > $out
	        done

		#Run intersection betweeen ATAC-seq files and all replicates
		bedtools multiinter -i *_sorted.bed ../../../ATAC-seq/dexamethasone/dexamethasone_consensuspeaks_sorted.bed -header > $output_dex
		bedtools multiinter -i *_sorted.bed ../../../ATAC-seq/control/control_consensuspeaks_sorted.bed -header > $output_control


		#Get corresponding gene for each peak
		output_dex_named="${output_dex%%.bed}_named.bed"
		bedtools intersect -a ../../../../TF_TSS_10kb_named.bed -b $output_dex > $output_dex_named
		output_control_named="${output_control%%.bed}_named.bed"
                bedtools intersect -a ../../../../TF_TSS_10kb_named.bed -b $output_control > $output_control_named

		#Cleanup
		rm *sorted.bed
		gzip ENC*.bed
  		cd ../
	done
	cd ../
done
