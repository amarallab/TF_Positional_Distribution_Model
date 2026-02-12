cd TF_Positional_Distribution_Model/data/ENCODE/RNAPol_Occupancy/

#Load deeptools
module load deeptools

#Compute the POLII occupancy  signal intensity for biologically relevant and spurious target genes
computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R True_TSS_TES.bed Spurious_TSS_TES.bed -S ENCFF155EYW.bigWig  -o output_matrix.gz
computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R True_sampled.bed Spurious_sampled.bed Missing_le2kb.bed -S ENCFF155EYW.bigWig  -o output_matrix_2kb.gz



#Plot signal profile
plotProfile -m output_matrix_2kb.gz --colors green purple blue -out ../../../Figures/Figure3B.png

#Plot heatmap
plotHeatmap -m output_matrix.gz \
-out ../../../Figures/Figure3Ctop.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -4 --zMax 4  

plotHeatmap -m output_matrix_2kb.gz \
-out ../../../Figures/Figure3Cbottom.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -4 --zMax 4  

