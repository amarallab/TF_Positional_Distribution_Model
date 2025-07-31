#Load required modules
import os 
import pandas as pd
import shutil
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

#Set working directory
working_dir ="TF_Positional_Distribution_Model/data/ENCODE/RNA-seq"
os.chdir(working_dir)

#Get list of gene IDs
id_all = data.gene_ids()

#Run compilation for all cell lines
compiled = pd.DataFrame()
compiled['Gene'] = data.gene_names()
files = os.listdir()
files = [f for f in files if '.tsv' in f]


for file in files:

    #Check file format and load accordingly
    with open(file) as f:
        first_line = f.readline()
    if '#' in first_line:
        df = pd.read_table(file,  skiprows=1)
        df['gene_id'] = df.Geneid
        df['RPK'] = df[df.columns[-3]]/df.Length
        df['TPM'] = df['RPK']/(df.RPK.sum()*1000000)
    else:
        df = pd.read_table(file)

    #Extract gene name from ENSEMBL gene ID
    df['gene_id'] = [i.split('.')[0] for i in df['gene_id']]
    df = df.loc[df.gene_id.isin(id_all)]
    df['Gene'] = [data.gene_by_id(i).gene_name for i in df.gene_id]

    #Extract TPM and compile with all samples in cell line
    df = df[['Gene','TPM']]
    df = df.groupby('Gene')['TPM'].mean().reset_index(name = 'TPM')
    compiled = compiled.merge(df, on='Gene')
    
compiled.columns = ['Gene'] + ['K562_' + str(i)  for i in range(1,len(compiled.columns))]
compiled.to_csv('Compiled_RNAseq.csv')
