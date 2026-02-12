import os
import pandas as pd
import tables
import numpy as np

training_K562 = pd.read_hdf('training.h5',
'training').set_index(['enhancer_name', 'promoter_name'])
training_K562 = training_K562.loc[training_K562.label == 1]

promoter_bed = training_K562.reset_index()[['promoter_chrom','promoter_start', 'promoter_end']].drop_duplicates()
promoter_bed.to_csv('K562_EnhancerAssociatedPromoters.bed',
                   header = None, index = None, sep = '\t')

enhancer_bed = training_K562.reset_index()[['enhancer_chrom','enhancer_start', 'enhancer_end']].drop_duplicates()
enhancer_bed.to_csv('K562_PromoterAssociatedEnhancer.bed',
                   header = None, index = None, sep = '\t')