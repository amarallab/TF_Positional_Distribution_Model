#Load required modules
import os 
import pandas as pd
import shutil
from pyensembl import EnsemblRelease

# release 77 uses human reference genome GRCh38
data = EnsemblRelease(77)

#Set working directory
working_dir ="/projects/b1042/AmaralLab/Maalavika/TF_Positional_Distribution_Model/data/ENCODE/A549_Perturbation/RNA-seq"
os.chdir(working_dir)

#Load metadata and extract relevant information
expt = pd.read_table("metadata.tsv")
# expt = expt.loc[(expt['Output type'] == 'gene quantifications') & (expt['File assembly'] == 'GRCh38')&
#                (expt['Lab'] == 'ENCODE Processing Pipeline')]

expt = expt[['Biosample term name', 'Assay', 'File accession']]


#Read each file (accessionID.tsv) and store counts and tpm values under cellline label
counts = pd.DataFrame()
tpm = pd.DataFrame()
colnames =[]
count = []

for treatment in ['control','dexamethasone']:
    os.chdir(treatment)
    
    for idx,i in enumerate(expt['File accession']):
        encode_files = os.listdir()
        encode_files = [i.replace('.tsv','') for i in encode_files if '.tsv' in i]
        
        if i in encode_files:
            df =pd.read_csv(i+'.tsv', index_col = None, sep = '\t',comment='#')
            df['counts'] = df[df.columns[-1]]
            df['gene_id'] = df.Geneid
            df['RPK'] = df['counts']/df.Length
            df['TPM'] = df['RPK']/(df.RPK.sum()*1000000)
            
            #Extract gene name from ENSEMBL gene ID
            df['gene_id'] = [i.split('.')[0] for i in df['gene_id']]
            df = df.loc[df.gene_id.isin(id_all)]
            df['Gene'] = [data.gene_by_id(i).gene_name for i in df.gene_id]

            #Extract TPM and compile with all samples in cell line
            tpm_df = df[['Gene','TPM']].groupby('Gene')['TPM'].mean().reset_index(name = 'TPM')
            tpm_df.index = tpm_df.Gene
            tpm = pd.concat([tpm, tpm_df['TPM']], axis = 1)
            counts_df = df[['Gene','counts']].groupby('Gene')['counts'].mean().reset_index(name = 'counts')
            counts_df.index = counts_df.Gene
            counts= pd.concat([counts, counts_df['counts']], axis = 1)
            count.append(colnames.count(treatment)+1)
            colnames.append(treatment)

    os.chdir('../')

#Clean up data - drop genes with missing values
tpm = tpm.dropna()
counts = counts.dropna()

#rename samples based on cell line
names = [str(i)+'_'+str(j) for i,j in zip(colnames, count)]
counts.columns = names
tpm.columns = names

#Save counts and TPM matrices
os.mkdir('processed')
counts.to_csv('processed/A549_dex12h_RNAseq_counts.csv')
tpm.to_csv('processed/A549_dex12h_RNAseq_TPM.csv')




