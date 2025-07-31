#!/bin/bash
#SBATCH -A b1042 ## <-- EDIT THIS TO BE YOUR ALLOCATION
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --partition  genomics
#SBATCH --job-name=targetGene
#SBATCH --mail-type=END,FAIL

module load bedtools #peak calling

#Set working directory
cd TF_Positional_Distribution_Model/data/Goode_ESC/ 

#Run each bed file in a given cell line to get target peaks

#For mm10 reference genome alignment
cd mm10
for i in $(ls -d */); do
	cd $i
	for bed_file in $(ls *narrowPeak); do
		save="${bed_file%peaks_peaks.narrowPeak}_targets_cds.bed"
       	        bedtools intersect -a ../TF_TSS_cds_named_mm10.bed -b $bed_file > $save

		save="${bed_file%peaks_peaks.narrowPeak}_targets_10kbp.bed"
                bedtools intersect -a ../TF_TSS_10kbp_named_mm10.bed -b $bed_file > $save

	done
	cd ../
done


cd ../mm39
for bed_file in $(ls *narrowPeak); do
	save="${bed_file%peaks_peaks.narrowPeak}_targets_cds.bed"
        bedtools intersect -a TF_TSS_cds_named_mm39.bed -b $bed_file > $save
        
        save="${bed_file%peaks_peaks.narrowPeak}_targets_10kbp.bed"
        bedtools intersect -a TF_TSS_10kbp_named_mm39.bed -b $bed_file > $save
        
done  

