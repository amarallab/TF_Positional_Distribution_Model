This folder contains scripts to reproduce Figure 6. 


It uses the processed files for K562 to identify promoter and enhancer bond sites in the genome. Prior to running the code,create a folder under the TF_Positional_Distribution_Model/data/ENCODE/K562 called EPI.

First, download https://github.com/shwhalen/targetfinder/tree/master/paper/targetfinder/K562/output-ep/training.h5 in this folder.


Next, run intersect.sh to get TF binding sites within annotated enhancer and promoter regions. 

Run Figur6.ipynb to generate plots for figure 6. 
