Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

Repository associated to the publication:

N-glycan signatures identified in tumor interstitial fluid and serum of breast cancer patients - association with tumor biology and clinical outcome
Thilde Terkelsen, Vilde D Haakensen, Radka Saldova, Pavel Gromov, Merete Kjær Hansen, Henning Stöckmann, Ole Christian Lingjærde, Anne-Lise Børresen-Dale, Elena Papaleo, Åslaug Helland, Pauline M. Rudd, Irina Gromova
submitted.


This repository contains N-glycan abundance data from tumour/normal interstitial fluids and serum obtained from a cohort of patients with breast cancer. The repository was made with intent of openly sharing both data and R-scripts used in analysis in relation to the publication.

The repository contains two folders:
                                    
          (1) N-glycan abundance data and patient metadata. 
          (2) Rscripts:
                           (a) N_Glycan_Functions.R
                           (b) N_Glycan_Analysis_All.R
                           (c) N_Glycan_Analysis_Paired.R
                           (d) N_Glycan_Survival_Analysis.R
                           (e) N_Glycan_Serum_Analysis.R
                           (f) N_Glycan_Figures.R
                                    
Aditionally the folder with R-scrips contains a sub-folder with dataframes geneated by scripts a-e, which may be used for easy plotting using script f, or to check results.

Requirements:
R version 3.3.1 or high
Rstudio version 1.1.383 or higher

Although R-packages should automatically be installed and errors raised if they cannot be, we here provide the user with the list of required packages:

CRAN:
limma
sva
openxlsx
ggplot2
dendextend
heatmap.plus
reshape
gdata
plyr
data.table
RColorBrewer
squash
survminer
car

Bioconductor:
survcomp


Notes:
a) We suggest to use Rstudio to run the scripts of interest. 
b) The scripts need to be run in the same folder where the N-glycan data and patient metadata are stored. Otherwise, the user would need to modify them with the right path.
b) N:B it is vital that the script N_Glycan_Functions.R is always run as the initital script as this scripts contains packages and costum functions needed for running the rest of the code.
