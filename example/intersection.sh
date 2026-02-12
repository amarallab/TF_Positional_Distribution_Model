module load bedtools

file='example.bed' ##Add path to bed file
ref_genome='reference_genome.bed' ##Add path to reference genome

output="${file%%.bed}_named.bed"

bedtools intersect -a $ref_genome -b $file > $output
