
# N-glycan signatures identified in tumor interstitial fluid and serum of breast cancer patients

Date: 06-02-2018
Author of R-scrips: Thilde Bagger Terkelsen, PhD student CBL, DCRC
Email: thilde@cancer.dk 


This repository contains N-glycan abundance data from tumour/normal interstitial fluids and serum obtained from a cohort of patients with breast cancer. The repository was made with intent of openly sharing both data and R-scripts used in analysis in relation to the publication "N-glycan signatures identified in tumor interstitial fluid and serum of breast cancer patients".

The repository cotains two folders:
                                    
          (1) N-glycan abundance data and patient metadata. 
          (2) Rscripts:
                           (a) N_Glycan_Functions.R
                           (b) N_Glycan_Analysis_All.R
                           (c) N_Glycan_Analysis_Paired.R
                           (d) N_Glycan_Survival_Analysis.R
                           (e) N_Glycan_Serum_Analysis.R
                           (f) N_Glycan_Figures.R
                                    
Aditionally the folder with R-scrips contins a sub-folder with dataframes geneated by scripts a-e, which may be used for easy plotting using script f, or to check results.

It is vital that the script N_Glycan_Functions.R is always run as the initital script as this scripts contains packages and costum functions needed for running the rest of the code.
