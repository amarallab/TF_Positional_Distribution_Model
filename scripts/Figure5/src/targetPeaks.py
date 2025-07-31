import os
import pandas as pd
import numpy as np



#Set working directory
working_dir ="TF_Positional_Distribution_Model/data/ENCODE/A549_Perturbations/ATAC-seq"
os.chdir(working_dir)

#Load files
dex_df = pd.read_table('dexamethasone_12h/dexamethasone_12h_consensuspeaks.bed')
control_df = pd.read_table('control/control_consensuspeaks.bed')

#select only peaks overlapping with target region
dex_df = dex_df.loc[(dex_df.num>1) & (dex_df['../../../../TF_TSS_10kb_named.bed'] ==1),dex_df.columns[:-1]]
control_df =control_df.loc[(control_df.num > 1) & (control_df['../../../../TF_TSS_10kb_named_updated.bed'] ==1), control_df.columns[:-1]]

#Estimate concordance
dex_df['concordance'] = (dex_df.num -1)/(len(dex_df.columns)-5)
control_df['concordance'] = (control_df.num -1)/(len(control_df.columns)-5)

os.mkdir('processed')

#Save
dex_df.to_csv('processed/dexamethasone_12h_consensuspeaks.bed', sep = '\t', index = False)
control_df.to_csv('processed/control_consensuspeaks.bed', sep = '\t', index = False)
dex_df.loc[dex_df.concordance ==1].to_csv('processed/dexamethasone_12h_consensuspeaks_concordance_1.bed', sep = '\t', index = False, header = False)
control_df.loc[control_df.concordance ==1].to_csv('processed/control_consensuspeaks_concordance_1.bed', sep = '\t', index = False, header = False)



