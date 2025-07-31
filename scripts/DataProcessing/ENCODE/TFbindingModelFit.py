#Load required modules
import os 
import pandas as pd
import numpy as np
from plotnine import *
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr, fisher_exact
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from math import exp,log
from scipy.stats import norm
from utils import *
from sklearn.metrics import r2_score
import patchworklib as pw
from random import choices
import math


#Set working directory
working_dir ="TF_Positional_Distribution_Model/data/ENCODE/"
os.chdir(working_dir)

#15 TFs used for the analysis
tfs = ['CEBPB', 'CTCF', 'ELK1', 'EP300', 'ESRRA', 'HDAC2', 'MAFK', 'MAX', 'MAZ', 'MYC', 'RAD21', 'RCOR1', 'REST', 'RFX5', 'SIN3A']
cellLines = ['A549','GM12878','HepG2','K562','MCF-7']

#Initialize lists/dataframes
params =[]
stats= []
compiled_cdf = pd.DataFrame()


#Get reference genome 
TSS_ref = pd.read_table('TF_TSS_10kb_named.bed')
TSS_ref['TSS'] = TSS_ref['end'] - 10000

for cl in cellLines:
    for idx,tf in enumerate(tfs):
    
        #Load file for TF in cell linecat TF
        df= pd.read_table(cl +'/ChIP-seq/'+tf+ '_consensus_accessible_peaks.bed')
        named_df =  pd.read_table(cl +'/ChIP-seq/'+tf+ '_consensus_accessible_peaks_named.bed',
                                 names = ['chrom','start','end','Gene'])

        #Merge to get named peaks
        df = df.merge(named_df, on=['chrom','start','end'])

        #Get distance of peak from TSS
        df=df.merge(TSS_ref[['TSS','Gene']])
        df['start_diff'] =  df['start'] - df['TSS']
        df['end_diff'] =  df['end'] - df['TSS']
        df['min_dist'] = df.apply(lambda x: abs(sum([x['start_diff'], x['end_diff']])/2), axis =1)

        #Select peaks that are accessible and present in atleast one replicate
        sel_column =[i for i in df.columns if tf in i]
        df = df.loc[(df[sel_column].sum(axis =1) > 0) & (df['../ATAC-seq/consensuspeaks_all_targets_10kb_sorted.bed']==1)]

        
        #Get cumulative count
        summary = [[i, sum(df.min_dist <= i)] for i in range(0,10000,100)]
        summary = pd.DataFrame(summary, columns =['interval','cumulative'])
        xdata = summary.interval.values
        y = summary.cumulative.values
        
        popt1, pcov1 = curve_fit(model, xdata, y, p0 =[100,2000,0], maxfev = 4000)
        params.append(np.append(popt1, [tf,cl]))
        stats.append(np.append(pcov1, [tf,cl]))

        #Compile data
        summary['TF'] = tf
        summary['CellLine'] =cl
        summary['Prediction'] = model(xdata, popt1[0], popt1[1],popt1[2])
        
        #Save fits (Optional)
##        p_log10 = (ggplot(summary)
##        +geom_point(aes(x = 'interval', y = 'cumulative'), alpha = 0.2)
##        +geom_line(aes(x='interval', y = 'Prediction'))
##        +ggtitle(tf + ' in ' + cl + ' - Exponential decay + constant spurious binding')
##        +theme_minimal()
##        +scale_x_log10())
##        p = (ggplot(summary)
##        +geom_point(aes(x = 'interval', y = 'cumulative'), alpha = 0.2)
##        +geom_line(aes(x='interval', y = 'Prediction'))
##        +ggtitle(tf + ' in ' + cl+ ' - Exponential decay + constant spurious binding')
##        +theme_classic())
##        filename =tf + '_' + cl+'.png'
##        ggsave(plot = p_log10, filename = filename, path = "Figures/log10_x/")
##        ggsave(plot = p, filename = filename, path = "Figures/original/")
        
        
        #Compile data
        compiled_cdf = pd.concat([compiled_cdf, summary])
      

#Create directory to store processed data
os.mkdir('processed')
os.chdir('processed')
#Save cumulative positional distribution of TFs and parameters for model fit 
compiled_cdf.to_csv('Compiled_accessible_target_ChIP_peaks.csv')
params = pd.DataFrame(params, columns = ['rb','tau','rs','TF','CellLine']).to_csv('ExpDecayConst_params.csv')

