

module load bedtools

#Get genome-wide accesible TF binding sites for identifying TF-bound enhancers
cd TF_Positional_Distribution_Model/data/ENCODE/K562

#Sort ATAC-seq file for multiinter (use file with all accessible sites, not just those within 10kbp of gene)
sort -k 1,1 -k2,2n ATAC-seq/consensuspeaks_all.bed > ATAC-seq/consensuspeaks_all_sorted.bed
	
#For each TF get set of consensus accessible peaks
cd ChIP-seq
for tf in 'CEBPB' 'CTCF' 'ELK1' 'EP300' 'ESRRA' 'HDAC2' 'MAFK' 'MAX' 'MAZ' 'MYC' 'RAD21' 'RCOR1' 'REST' 'RFX5' 'SIN3A'; do

	#Unzip each raw bed file
	gzip -d ENC*$tf*

	#Sort file
	for bed in $(ls ENC*$tf*); do
		out="${bed%%.bed}_sorted.bed"
		sort -k 1,1 -k2,2n $bed  > $out
	done

	#Run multiinter
	output="${tf}_consensus_accessible_peaks_genomewide.bed"
	bedtools multiinter -i *$tf*sorted.bed ../ATAC-seq/consensuspeaks_all_sorted.bed -header > $output
		
		
	#Remove intermediate sort files
	rm *sorted.bed
	gzip ENC*$tf*

done


#Create folder to save intersections
cd ../EPI
mkdir intersections

#Extract enhancer and promoter bed files
python ../../../../scripts/Figure6/src/ep_bedfiles.py

#Get genes within 10kb of annotated promoters
bedtools intersect -a TF_TSS_10kb.bed -b K562_EnhancerAssociatedPromoters.bed -wb -wa > K562_EnhancerAssociatedPromoters_named.bed 


#Run intersection for accessible peaks
for i in ../ChIP-seq/*_consensus_accessible_peaks.bed; do 

	filename="${i##*/}"
	output="intersections/${filename%%.bed}_intersection_promoters.bed"
	bedtools intersect -a $i -b K562_EnhancerAssociatedPromoters.bed -wa -wb > $output

	input="${i%%.bed}_genomewide.bed"
	output="intersections/${filename%%.bed}_intersection_enhancers.bed"
        bedtools intersect -a $input -b K562_PromoterAssociatedEnhancer.bed -wa -wb > $output


done

