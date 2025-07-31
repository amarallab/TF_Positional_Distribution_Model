
#Using release 69, same as GenBank GCA_000001635.9 (mm10)
#GCA_000001405.15 release 77 homo sapiens



#Load required modules
import os
import pandas as pd
from pyensembl import EnsemblRelease
import numpy as np 


#Set working directory
working_dir ="TF_Positional_Distribution_Model/data/ENCODE/"
os.chdir(working_dir)


#List functional chromosomes, used for analyses
chrom = ['chr'+str(i) for i in range(1,23)]
chrom = chrom + ['chrX', 'chrY']

#Create bed file for all TSS start sites
data = EnsemblRelease(77, species = 'human')
genes = data.gene_names()
gene_info = [['chr'+data.loci_of_gene_names(gene)[0].contig, #chromosome
              data.loci_of_gene_names(gene)[0].start ,  #Start of cds
              data.loci_of_gene_names(gene)[0].end,     #End of cds
              data.loci_of_gene_names(gene)[0].strand,  #Strand
              gene ] for gene in genes if gene != '']   #Name of gene
gene_info = pd.DataFrame(gene_info)

#Keep genes present on functional chromosomes
keepseq = [i in chrom for i in gene_info[0]]
gene_info = gene_info.loc[keepseq]

#Set TSS for genes on negative strand and end of cds
gene_info.columns = ['chrom','start','end','strand','Gene']
gene_info['TSS'] = [i.end if i['strand'] == '-' else i.start for idx, i in gene_info.iterrows()]


#Create reference annotation file with 10kbp +/- TSS
gene_info['start']= gene_info.TSS-10000
gene_info['end']= gene_info.TSS+10000

#Save reference annotation file
gene_info[['chrom','start','end','Gene']].to_csv('TF_TSS_10kb_named.bed', sep="\t",index=None )
