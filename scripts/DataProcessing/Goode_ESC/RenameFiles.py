#Load required modules
import os
import pandas as pd

#Set working director
working_dir="TF_Positional_Distribution_Model/data/Goode_ESC/"
os.chdir(working_dir)

#Extract metadata files for input and ChIP-seq data
metadata = os.listdir()
metadata = [f for f in metadata if 'RunTable.csv' in f]


#Rename cell types
conversion ={'Mesoderm': 'MES', 'Hemangioblast': 'HB', 'Hemogenic Endothelium':'HE', 'Hematopoietic Progenitors':'HP',
            'Macrophages': 'MAC'}

#Get new file names
files_old = os.listdir()
files_old = [f for f in files_old if 'bam' in f]
for f in metadata:
    meta_table = pd.read_csv(f)
    meta_table.loc[meta_table['ChIP_antibody'].astype('str') == 'nan', 'ChIP_antibody'] = 'input'
    meta_table['ChIP_antibody'] =[f.split(' ')[0] for f in meta_table['ChIP_antibody']]
    meta_table['Label'] = meta_table['source_name'].replace(conversion)
    meta_table['Label']=meta_table['Label'] +'_'+ meta_table['ChIP_antibody']
    g = meta_table.groupby(['Label'])
    meta_table.loc[g['Label'].transform('size').gt(1),
           'Label'] += '_'+ (g.cumcount()+1).astype(str)
    for idx,row in meta_table[['Run', 'Label']].iterrows():
        files_new = [i.replace(row[0],row[1]) for i in files_old] 
        files_old = files_new

#Rename all files 
files_old = os.listdir()
files_old = [f for f in files_old if 'bam' in f]
for idx,file in enumerate(files_old):
    os.rename(file, files_new[idx])
