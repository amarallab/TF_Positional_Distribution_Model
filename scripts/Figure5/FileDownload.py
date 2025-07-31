#Load required modules
import os
import pandas as pd
import numpy as np
import shutil
import urllib.request


#Set working directory
working_dir ="TF_Positional_Distribution_Model/data/ENCODE/A549_Perturbation"
os.chdir(working_dir)

## Downlaod ChIP-seq data
os.chdir('ChIP-seq')
df = pd.read_table('metadata.tsv')

#Select IDR thresholded peaks and dexamethasone treatment samples
df_sub = df.loc[(df['Output type'] == 'IDR thresholded peaks') &
               (df['Biosample term name'] == 'A549') &
               ((df['Biosample treatments duration'] == '12 hour') |
                (df['Biosample treatments'].isnull())) &
                ([',' not in str(i) for i in df['Biological replicate(s)']])]
df_sub['Biosample treatments'] = df_sub['Biosample treatments'].replace(np.nan, 'control')

# Get set of TFs present in control and dexamethasone treatment
tfs = df_sub.groupby(['Experiment target'])['Biosample treatments'].apply(lambda x: len(np.unique(x)) >1 )
tfs= tfs.index[tfs]

# Subset files with TFs
df_sub =df_sub.loc[df_sub['Experiment target'].isin(tfs)]

#Download files
for idx, line in df_sub.reset_index().iterrows():
    tf = line['Experiment target'].split('-')[0]
    path = line['Biosample treatments']
    if not os.path.exists(path):
        os.mkdir(path)
        
    path = path + '/' +tf
    if not os.path.exists(path):
        os.mkdir(path)
    file = 'https://www.encodeproject.org/files/' + line['File accession']+ '/@@download/' + line['File accession']+ '.bed.gz'
    savefile = path  + '/' + line['File accession']  + '.bed.gz'
    urllib.request.urlretrieve(file, savefile)



## Download ATAC-seq data
os.chdir('../ATAC-seq')
df = pd.read_table('metadata.tsv')


#Select files corresponding to IDR thresholded peaks, and dex12h samples
df_sub = df.loc[(df['Output type'] == 'IDR thresholded peaks') &
               (df['Biosample term name'] == 'A549') &
               ((df['Biosample treatments duration'] == '12 hour') |
               (df['Biosample treatments'].isnull())) &
                ([',' not in str(i) for i in df['Biological replicate(s)']])]
df_sub['Biosample treatments'] = df_sub['Biosample treatments'].replace(np.nan, 'control')

#Download files
for idx, line in df_sub.reset_index().iterrows():
    path = line['Biosample treatments']
    if not os.path.exists(path):
        os.mkdir(path)
        
    file = 'https://www.encodeproject.org/files/' + line['File accession']+ '/@@download/' + line['File accession']+ '.bed.gz'
    savefile = path  + '/' + line['File accession']  + '.bed.gz'
    urllib.request.urlretrieve(file, savefile)



## Download RNA-seq data
os.chdir('../RNA-seq')
df = pd.read_table('metadata.tsv')


df['Biosample treatments'] = df['Biosample treatments'].replace(np.nan, 'control')

#Download files
for idx, line in df.reset_index().iterrows():
    path = line['Biosample treatments']
    if not os.path.exists(path):
        os.mkdir(path)
        
    file = 'https://www.encodeproject.org/files/' + line['File accession']+ '/@@download/' + line['File accession']+ '.tsv'
    savefile = path  + '/' + line['File accession']  + '.tsv'
    urllib.request.urlretrieve(file, savefile)

