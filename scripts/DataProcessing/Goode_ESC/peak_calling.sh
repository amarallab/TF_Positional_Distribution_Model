#!/bin/bash
#SBATCH -A b1042 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --partition  genomics
#SBATCH --job-name=peak-calling
#SBATCH --mail-type=END,FAIL
###Step5 peak-calling

module load MACS2/2.1.0 #peak calling

cd TF_Positional_Distribution_Model/data/Goode_ESC/mm10

#Run peak calling fro multiple p-values
for i in $(ls *nodup.bam | sed -e 's/\_.*$//' |sort -u);
do
  	input="${i}_input.sorted.nodup.bam"
        for j in $(ls $i*nodup.bam);
        do
          	for ((val=1; val<10; val++));
                do
                  	pval=$(echo "10^-($val)"|bc -l)
                        folder="1e-$val"
                        
                        #Create directory for each p-value threshold
                        mkdir $folder
                        j_out=${j%".sorted.nodup.bam"}
                        output="${folder}/"$j_out"_1e-$val"
        
                        #Call peaks
                        macs2 callpeak -t $j  -c $input\
                        -f BAM -n $output -g hs -p $pval  --nomodel
                        echo $j completed

                done
        done
done            

cd ../mm39
for i in $(ls *nodup.bam | sed -e 's/\_.*$//' |sort -u);
do
  	input="${i}_input.sorted.nodup.bam"
        for j in $(ls $i*nodup.bam);do
                j_out=${j%".sorted.nodup.bam"}
                        
                #Call peaks
                macs2 callpeak -t $j  -c $input\
                -f BAM -n $j_out -g hs -p 1e-5  --nomodel
                echo $j completed
                
               	done
        done
done



