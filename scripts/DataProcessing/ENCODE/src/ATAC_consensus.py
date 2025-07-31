#Load required modules
import os
import pandas as pd
import numpy as np

#Load files
df = pd.read_table('consensuspeaks.bed')

#Estimate concordance
df['concordance'] = (df.num)/(len(df.columns)-5)

#Save peaks present in all replicates
df.loc[df.concordance ==1].to_csv('consensuspeaks_all.bed', sep = '\t', index = False, header = False)



