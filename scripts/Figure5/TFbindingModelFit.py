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
working_dir ="TF_Positional_Distribution_Model/data/ENCODE/A549_Perturbation/ChIP-seq"
os.chdir(working_dir)

#15 TFs used for the analysis
treatments = ['control','dexamethasone']

#Initialize lists/dataframes
params =[]
stats= []
compiled_cdf = pd.DataFrame()

#Get reference genome 
TSS_ref = pd.read_table('../../TF_TSS_10kb_named.bed')
TSS_ref['TSS'] = TSS_ref['end'] - 10000


for idx,treatment in enumerate(treatments):

    #List all TFs in each condition
    tfs= os.listdir(treatment)

    #Run fits for all TFs
    for tf in tfs:
    
        #Load file for TF in cell line

        df =  pd.read_table(treatment +'/'+tf + '/'+ treatment + '_consensuspeaks_accessible_'+ treatment +'.bed')
        named_df =  pd.read_table(treatment +'/'+tf + '/'+ treatment + '_consensuspeaks_accessible_'+ treatment +'_named.bed',
                                     names = ['chrom','start','end','Gene' ])

        #Merge to get named peaks
        df = df.merge(named_df, on=['chrom','start','end'])

        #Get distance of peak from TSS
        df=df.merge(TSS_ref[['TSS','Gene']])
        df['start_diff'] =  df['start'] - df['TSS']
        df['end_diff'] =  df['end'] - df['TSS']
        df['min_dist'] = df.apply(lambda x: abs(sum([x['start_diff'], x['end_diff']])/2), axis =1)

        #Select peaks that are accessible and present in atleast one replicate
        sel_column =[i for i in df.columns if 'ENC' in i]
        atac_col = [i for i in df.columns if 'ATAC-seq' in i][0]
        df = df.loc[(df[sel_column].sum(axis =1) > 0) & (df[atac_col]==1)]

        
        #Get cumulative count
        summary = [[i, sum(df.min_dist <= i)] for i in range(0,10000,100)]
        summary = pd.DataFrame(summary, columns =['interval','cumulative'])
        xdata = summary.interval.values
        y = summary.cumulative.values
        
        popt1, pcov1 = curve_fit(model, xdata, y, p0 =[100,2000,0], maxfev = 4000)
        params.append(np.append(popt1, [tf,treatment]))
        stats.append(np.append(pcov1, [tf,treatment]))

        #Compile data
        summary['TF'] = tf
        summary['Treatment'] = treatment
        summary['Prediction'] = model(xdata, popt1[0], popt1[1],popt1[2])
        

        #Compile data
        compiled_cdf = pd.concat([compiled_cdf, summary])
      

#Create directory to store processed data
os.mkdir('processed')
os.chdir('processed')
#Save cumulative positional distribution of TFs and parameters for model fit 
compiled_cdf.to_csv('Compiled_accessible_target_ChIP_peaks.csv')
params = pd.DataFrame(params, columns = ['rb','tau','rs','TF','Treatment']).to_csv('ExpDecayConst_params.csv')

