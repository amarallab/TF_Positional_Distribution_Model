This folder contains the scripts needed to generate Figure 3 of the paper.

## Downloading the data
To download the data from ENCODE, run the following line of code in the 
TF_Positional_Distribution_Model/data/ENCODE/RNAPol_Occupancy/ directory

wget https://www.encodeproject.org/files/ENCFF155EYW/@@download/ENCFF155EYW.bigWig

This downloads the RNA Polymerase 2 subunit A (POLR2A) occupancy data in .bigWig format.

Output: ENCFF155EYW.bigWig

## Classifying genes as biologically relevant and spurious targets of MYC (RNApol2_tracks.py)
This script classifies binding events for MYC as biologically relevant and spurious. For each class, a reference annotation
file containing the target gene loci (chromosome, TSS, TES, strand) and the name is created and saved as a .bed file. 
 
Output: True_TSS_TES_MYC.bed, Spurious_TSS_TES_MYC.bed 


## Generating Figure 3 (Figure3.sh)
This bash script generates the panels in Figure 3 in the paper using deeptools. A heatmap and a signal intensity profile 
for POLII occupancy at the TSS of biologically relevant and spurious target genes are generated. 

Output: TF_Positional_Distribution_Model/Figures/Figure3B-C.png


