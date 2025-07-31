#!/bin/bash
#SBATCH -A b1042 ## <-- EDIT THIS TO BE YOUR ALLOCATION
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --partition  genomics
#SBATCH --job-name=download
#SBATCH --mail-type=END,FAIL

#In this code we select IDR thresholded peak files for all ChIP-seq datasets from ENCODE that correspond to unperturbed cells


cd TF_Positional_Distribution_Model/scripts/DataProcessing/ENCODE/src
python fileSelection.py


