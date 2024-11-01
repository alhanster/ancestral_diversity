# Increased representation of human genetic ancestries improves genic intolerance metrics 


Code repository for calculating RVIS, MTR, LOF O/E, LOF-FDR scores.
These intolerance metrics were calculated using using two large whole-exome sequencing datasets: the UK Biobank (n=460,551) and gnomAD v2.1 (n=125,748). <br>
<br>
For questions, contact alh16@rice.edu.

# Usage
After downloading the repository, expand the file and run the following command lines. Tables and figures will be added to the output folder.
```
cd ancestry_diversity
Rscript figures/Figure1_RVIS.R
Rscript figures/Figure2_MTR.R
Rscript figures/Figure4_LOF.R
```

# Dependencies
This package depends on the following R packages:

install.packages("ggplot2") <br>
install.packages("patchwork") <br>
install.packages("tidyverse") <br>
install.packages("data.table") <br>
install.packages("R.utils") <br>
install.packages("pROC") <br>

Input: CSV files located in the data folder. The CSV files include variables necessary for calculating each intolerance metric.

# Figures

Figure1_RVIS.R computes ancestry specific RVIS in gnomAD and UKB. Logistic Regression is performed to evaluate how RVIS can be used to determine genes in disease implicated gene lists. DeLong test is used to determine whether difference between AUC values are statistically significant. <br>

Figure2_MTR.R computes MTR in UKB cohorts. Logistic Regression is performed to evaluate how MTR can be used to determine genes in disease implicated gene lists. DeLong test is used to determine whether difference between AUC values are statistically significant.  <br>

Figure4_LOF_OE.R computes LOF O/E in UKB cohorts. Logistic Regression is performed to evaluate how LOF O/E can be used to determine genes in disease implicated gene lists. DeLong test is used to determine whether difference between AUC values are statistically significant.  <br>

# Scripts

GeneList_VariantAnnotations.R includes variant annotations and pulls gene lists.

UKB_RVIS.R and gnomAD_RVIS.R include functions for calculating RVIS, logistic regression, and DeLong test. It also includes a function for graphing Figure 1 from our paper.

UKB_MTR.R and UKB_LOF.R include functions for logistic regression and DeLong test. It also includes a function for graphing Figure 2 and Figure 4 from our paper.

