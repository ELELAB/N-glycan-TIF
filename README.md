Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark.

Repository associated with the publication:

N-glycan signatures identified in tumor interstitial fluid and serum of breast cancer patients - association with tumor biology and clinical outcome
Thilde Terkelsen, Vilde D Haakensen, Radka Saldova, Pavel Gromov, Merete Kjær Hansen, Henning Stöckmann, Ole Christian Lingjærde, Anne-Lise Børresen-Dale, Elena Papaleo*, Åslaug Helland, Pauline M. Rudd, Irina Gromova, Mol Oncol, 2016, doi: 10.1002/1878-0261.12312

Corresponding author: Elena Papaleo, elenap@cancer.dk

R-scripts were created by: Thilde Bagger Terkelsen, thilde@cancer.dk

This repository contains N-glycan abundance data from tumor/normal interstitial fluids and serum obtained from a cohort of patients with breast cancer. The repository was made with intent of openly sharing both data and R-scripts used for analysis in relation to the publication.

The repository contains two folders:

    (1) N-glycan abundance data and patient metadata. These are the data used as the starting point for our analyses.
    (2) A collection of R-scripts that recapitulate our work:
                                    
         (a) N_Glycan_Functions.R (Packages and Custom Functions, N.B: run as the initial script)
         (b) N_Glycan_Analysis_All.R (Differential Expression Analysis using all TIF and NIF samples)
         (c) N_Glycan_Analysis_Paired.R (Differential Expression Analysis using only TIF-NIF pairs)
         (d) N_Glycan_Survival_Analysis.R (Survival Analysis, cox proportional hazard regression)
         (e) N_Glycan_Serum_Analysis.R (Serum-TIF Correlation Analysis, overlap with MDG cohort)
         (f) N_Glycan_Figures.R (Generates All Plots used in the publication)
                                    
Aditionally the folder with R-scripts contains a sub-folder with dataframes generated by scripts a-e, which may be used for easy plotting using script f, or to check results.


Requirements:
                                   
    R version 3.3.1 or higher
    Rstudio version 1.1.383 or higher        
    Bioconductor version 3.6 or higher	

Although R-packages should automatically be installed and errors raised if they cannot be, we here provide the user with the list of required packages for manual installation:

CRAN:

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
    scales
    statmod
    
Bioconductor:

    limma
    sva
    survcomp               

NOTES:

a) We suggest to use Rstudio to run the scripts of interest so that you can follow the analyses one line at the time and digest the results.

b) The scripts need to be run in the same folder where the N-glycan data and patient metadata are stored. Otherwise, the user would need to modify them specifying the right path to the input data.

c) N.B.  It is vital that the script N_Glycan_Functions.R is always run as the initial script as this scripts contains packages and costum functions needed for running the rest of the code.
