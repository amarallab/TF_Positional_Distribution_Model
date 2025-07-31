#Load required modules
import os
import pandas as pd
from pyensembl import EnsemblRelease
import numpy as np 
#Using release 69, same as GenBank GCA_000001635.9


#Set working directory
working_dir="TF_Positional_Distribution_Model/data/Goode_ESC/mm10"
os.chdir(working_dir)

#List of chromosomes to consider
chrom = ['chr'+str(i) for i in range(1,20)]
chrom = chrom + ['chrX', 'chrY']

#Reference genome (mm10)
data = EnsemblRelease(69, species = 'mouse')
genes = data.gene_names()
gene_info = [['chr'+data.loci_of_gene_names(gene)[0].contig,
              str(data.loci_of_gene_names(gene)[0].start) ,
              str(data.loci_of_gene_names(gene)[0].end), 
              str(data.loci_of_gene_names(gene)[0].strand), 
              gene ] for gene in genes if gene != '']


gene_info = pd.DataFrame(gene_info, 
                        columns=['chrom','TSS','TES','strand', 'gene'])
keepseq = [i in chrom for i in gene_info['chrom']]
gene_info = gene_info.loc[keepseq]
gene_info['TSS_corr'] = [i.TES if i.strand=="-" else i.TSS for idx,i in gene_info.iterrows()]
gene_info['start'] = gene_info.TSS_corr.astype(int) - 10000
gene_info['end'] = gene_info.TSS_corr.astype(int) + 10000
gene_info[['chrom','start', 'end','gene']].to_csv('TF_TSS_10kbp_named_mm10.bed', sep = "\t", index = None)
gene_info[['chrom','TSS', 'TES','gene']].to_csv('TF_TSS_cds_named_mm10.bed', sep = "\t", index = None)

# Reference genome (mm39)
os.chdir('../mm39')
data = EnsemblRelease(109, species = 'mouse')
genes = data.gene_names()
gene_info = [['chr'+data.loci_of_gene_names(gene)[0].contig,
              str(data.loci_of_gene_names(gene)[0].start) ,
              str(data.loci_of_gene_names(gene)[0].end),
              str(data.loci_of_gene_names(gene)[0].strand),
              gene ] for gene in genes if gene != '']


gene_info = pd.DataFrame(gene_info,
                        columns =['chrom','TSS','TES','strand', 'gene'])
keepseq = [i in chrom for i in gene_info['chrom']]
gene_info = gene_info.loc[keepseq]
gene_info['TSS_corr'] = [i.TES if i.strand=="-" else i.TSS for idx,i in gene_info.iterrows()]
gene_info['start'] = gene_info.TSS_corr.astype(int) - 10000
gene_info['end'] = gene_info.TSS_corr.astype(int) + 10000

gene_info[['chrom','start', 'end','gene']].to_csv('TF_TSS_10kbp_named_mm39.bed', sep = "\t", index = None)
gene_info[['chrom','TSS', 'TES','gene']].to_csv('TF_TSS_cds_named_mm39.bed', sep = "\t", index = None)


