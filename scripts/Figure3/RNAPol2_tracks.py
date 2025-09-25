#Load required modules
import os
import pandas as pd
import numpy as np
from pyensembl import EnsemblRelease
from utils import *

#Set working directory
working_dir ="/projects/b1042/AmaralLab/Maalavika/TF_Positional_Distribution_Model/data/ENCODE/RNAPol_Occupancy"
os.chdir(working_dir)

#Reference genome
chrom = ['chr'+str(i) for i in range(1,23)]
chrom = chrom + ['chrX', 'chrY']
data = EnsemblRelease(77, species = 'human')
genes = data.gene_names()
gene_info = [['chr'+data.loci_of_gene_names(gene)[0].contig,
              data.loci_of_gene_names(gene)[0].start ,
              data.loci_of_gene_names(gene)[0].end,
              gene, 
              data.loci_of_gene_names(gene)[0].start ,
              data.loci_of_gene_names(gene)[0].strand ] for gene in genes if gene != '']

gene_info = pd.DataFrame(gene_info)
keepseq = [i in chrom for i in gene_info[0]]
ref = gene_info.loc[keepseq]
ref.columns = ['chrom','TSS','TES','Gene','TSS_corr','strand']
ref['TSS_corr'] = [i.TES if i.strand =="-" else i.TSS for idx,i in ref.iterrows()]

#Get named peaks
df= pd.read_table('../A549/ChIP-seq/MYC_consensus_accessible_peaks.bed')
df_named = pd.read_table('../A549/ChIP-seq/MYC_consensus_accessible_peaks_named.bed', names = ['chrom','start','end','Gene'])
df = df.merge(df_named, on=['chrom','start','end'])

#Select accessible peaks
df = df.loc[df['../ATAC-seq/consensuspeaks_all_targets_10kb_sorted.bed'] ==1]
sel_columns = [i for i in df.columns if 'MYC' in i ]
df= df.loc[df[sel_columns].sum(axis=1) >0]



#Get distance of peaks from TSS
df = df.merge(ref[['TSS_corr','Gene']], on = 'Gene')
df['start_diff'] =  df['start'] - df['TSS_corr']
df['end_diff'] =  df['end'] - df['TSS_corr']
df['min_dist'] = df.apply(lambda x: abs(sum([x['start_diff'], x['end_diff']])/2), axis =1)

#Estimate parameters
params = pd.read_csv('../processed/ExpDecayConst_params.csv',index_col =0)
params = params.loc[(params.CellLine=="A549") &
                   (params.TF == "MYC")]

#Calculate threshold within which true binding to spurious binding ratio is 10. 
t = n_max_true(params.rb.values[0], params.tau.values[0], params.rs.values[0], [0.1]).Distance

df['Threshold'] = df.min_dist <= t.values[0]

#Classify true and  spurious edges
classified_gene_list = df.groupby('Gene')['Threshold'].any().reset_index(name='present')
true_edges = ref.loc[ref.Gene.isin(classified_gene_list.Gene.loc[classified_gene_list.present])]
false_edges = ref.loc[ref.Gene.isin(classified_gene_list.Gene.loc[~classified_gene_list.present])]

#Save reference files for true and spurious edges
true_edges.to_csv('True_TSS_TES.bed', sep = '\t', index = False, header = False)
false_edges.to_csv('Spurious_TSS_TES.bed', sep = '\t', index = False, header = False)
