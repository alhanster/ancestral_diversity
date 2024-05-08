# Increased representation of human genetic ancestries improves genic intolerance metrics 


Code repository for calculating RVIS, MTR, LOF O/E, LOF FDR scores. We need to explain more about the code!

# Usage
After downloading the repository, expand the file and paste the path into setwd() in R files in figures folder.

# Dependencies
This package depends on the following R packages:

install.packages("ggplot2") <br>
install.packages("patchwork") <br>
install.packages("tidyverse") <br>
install.packages("data.table") <br>
install.packages("R.utils") <br>
install.packages("pROC") <br>

# Figures

Figure1_RVIS.R computes ancestry specific RVIS in gnomAD and UKB. Logistic Regression is performed to evaluate how RVIS can be used to determine genes in disease implicated gene lists. DeLong test is used to determine whether difference between AUC values are statistically significant. <br>

Figure2_MTR.R computes MTR in UKB cohorts. Logistic Regression is performed to evaluate how MTR can be used to determine genes in disease implicated gene lists. DeLong test is used to determine whether difference between AUC values are statistically significant.  <br>

Figure4_LOF.R computes LOF O/E and LOF FDR in UKB cohorts. Logistic Regression is performed to evaluate how LOF O/E and LOF FDR can be used to determine genes in disease implicated gene lists. DeLong test is used to determine whether difference between AUC values are statistically significant.  <br>


