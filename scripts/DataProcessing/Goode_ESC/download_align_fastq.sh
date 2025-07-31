#!/bin/bash
#SBATCH -A b1042 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --partition  genomics
#SBATCH --job-name=fastq-download
#SBATCH --mail-type=END,FAIL

module load bowtie2 #alignment
module load samtools

#Download input files
web="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc="
while read -r line; 
do 
	out="${line}.fastq.gz"
	site=$web$line
	wget -O $out $site
	gzip -d $out

	#Label output files
	out="${line}.fastq"
	
	###Align against mm10
	bam="mm10/$line.bam"
	bowtie2 -q $out  -x mm10\
	|samtools view  -bS - > $bam	

	#sort
	sbam="mm10/$line.sorted.bam"
	samtools sort $bam  -o $sbam 

	#remove pcr duplicates
	sdbam="mm10/$line.sorted.nodup.bam"
	samtools rmdup -s $sbam $sdbam
	

	###Align against mm39
        bam="mm39/$line.bam"
        bowtie2 -q $out  -x mm39/GRCm39/GRCm39\   
        |samtools view  -bS - > $bam
        
        #sort 
        sbam="mm39/$line.sorted.bam"
        samtools sort $bam  -o $sbam
        
        #remove pcr duplicates
        sdbam="mm39/$line.sorted.nodup.bam"
        samtools rmdup -s $sbam $sdbam  

	
	#remove any input files (fastq)
	rm *fastq
done < INPUT_SRR_Acc_List.txt

#Download TF ChIP experiments

web="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc="
while read -r line;
do
	out="${line}.fastq.gz"
        site=$web$line
        wget -O $out $site
        gzip -d $out

        #Label output files
        out="${line}.fastq"

        ###Align against mm10
        bam="mm10/$line.bam"
        bowtie2 -q $out  -x mm10\
        |samtools view  -bS - > $bam

        #sort
        sbam="mm10/$line.sorted.bam"
        samtools sort $bam  -o $sbam

        #remove pcr duplicates
        sdbam="mm10/$line.sorted.nodup.bam"
        samtools rmdup -s $sbam $sdbam


        ###Align against mm39
        bam="mm39/$line.bam"
        bowtie2 -q $out  -x mm39/GRCm39/GRCm39\
        |samtools view  -bS - > $bam

        #sort
        sbam="mm39/$line.sorted.bam"
        samtools sort $bam  -o $sbam

        #remove pcr duplicates
        sdbam="mm39/$line.sorted.nodup.bam"
        samtools rmdup -s $sbam $sdbam


        #remove any input files (fastq)
        rm *fastq

done < SRR_Acc_List.txt

