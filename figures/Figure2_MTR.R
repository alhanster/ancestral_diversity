################################################################################
# Purpose: Figure 2: UKB MTR                                                   #
# Author: Alexander Han                                                        #
################################################################################

# Paste Path to Repository
setwd()

# Pull Variant Annotation and Gene Lists
source("scripts/GeneList_VariantAnnotations.R")

# Pull Functions for UKB MTR Computation, Logistic Regression, DeLong Test
source("~/Desktop/Scripts/UKB_MTR.R")

# Loading UKB Data
data <- fread("data/UKB_MTR.csv")

df <- data |>
  select(-MTR_Diverse43k, -MTR_NFE43k, -MTR_NFE440k, -MTR_FullDataset, -MTR_AFR, -MTR_ASJ, -MTR_EAS, -MTR_SAS, -MTR_NFE20k)

# Computing MTR Score
MTR <- df %>%
  mutate(
    `Maximally Diverse (n=43k)` = (Diverse_mis/(Diverse_mis + Diverse_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `NFE (n=43k)` = (Nfe_43k_mis/(Nfe_43k_mis + Nfe_43k_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `NFE (n=440k)` = (Nfe_440k_mis/(Nfe_440k_mis + Nfe_440k_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `Full Dataset (n=460k)` = (All_mis/(All_mis + All_syn)) / (possible_mis/(possible_mis + possible_syn)),
    
    `AFR` = (afr_mis/(afr_mis + afr_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `ASJ` = (asj_mis/(asj_mis + asj_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `EAS` = (eas_mis/(eas_mis + eas_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `SAS` = (sas_mis/(sas_mis + sas_syn)) / (possible_mis/(possible_mis + possible_syn)),
    `NFE (n=20k)` = (nfe_20k_mis/(nfe_20k_mis + nfe_20k_syn)) / (possible_mis/(possible_mis + possible_syn))
  ) %>% 
  select(Gene, `Maximally Diverse (n=43k)`, `NFE (n=43k)`, `NFE (n=440k)`, `Full Dataset (n=460k)`, `AFR`, `ASJ`, `EAS`, `SAS`, `NFE (n=20k)`)

# Logistic Regression and DeLong Test
dee_results <- PrintAUC(MTR, score_cols, dee_monoallelic)
dee_AUC <- dee_results$auc_df %>% 
  rename(`DEE Monoallelic\n (n=94)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dee_delong <- dee_results$delong_results %>% 
  mutate(Gene_List = "DEE Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

dd_results <- PrintAUC(MTR, score_cols, dd_monoallelic)
dd_AUC <- dd_results$auc_df %>% 
  rename(`DD Monoallelic\n (n=435)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dd_delong <- dd_results$delong_results %>% 
  mutate(Gene_List = "DD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

asd_results <- PrintAUC(MTR, score_cols, asd_monoallelic)
asd_AUC <- asd_results$auc_df %>% 
  rename(`ASD Monoallelic\n (n=190)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
asd_delong <- asd_results$delong_results %>% 
  mutate(Gene_List = "ASD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

mgi_results <- PrintAUC(MTR, score_cols, mgi_essential)
mgi_AUC <- mgi_results$auc_df %>% 
  rename(`Mouse Essential\n (n=2454)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
mgi_delong <- mgi_results$delong_results %>% 
  mutate(Gene_List = "Mouse Essential") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

HI_results <- PrintAUC(MTR, score_cols, clingen_HI)
HI_AUC <- HI_results$auc_df %>% 
  rename(`Haploinsufficient\n (n=360)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
HI_delong <- HI_results$delong_results %>% 
  mutate(Gene_List = "Haploinsufficient") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)


# Compiling AUC Values for Logistic Regression
MTR_AUC <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
MTR_AUC <- MTR_AUC %>% 
  rename("Ancestry" = score)
MTR_AUC

delongtest <- rbind(dee_delong, dd_delong, asd_delong, mgi_delong, HI_delong)
delongtest

write.csv(delongtest, "MTR_Score_DeLong_test.csv", row.names = FALSE)

# Logistic Regression
dee_AUC <- PrintLogRegResults(MTR, score_cols, dee_monoallelic, maf.threshold = 0.0005) %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

dd_AUC <- PrintLogRegResults(MTR, score_cols, dd_monoallelic, maf.threshold = 0.0005) %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

asd_AUC <- PrintLogRegResults(MTR, score_cols, asd_monoallelic, maf.threshold = 0.0005) %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

mgi_AUC <- PrintLogRegResults(MTR, score_cols, mgi_essential, maf.threshold = 0.0005) %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

HI_AUC <- PrintLogRegResults(MTR, score_cols, clingen_HI, maf.threshold = 0.0005) %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

# Compiling Logistic Regression
UKB_MTR_log <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
UKB_MTR_log <- UKB_MTR_log %>% 
  rename("Ancestry" = score)

# Compiling Figure
a <- PrintGraph(MTR_AUC, "DEE Monoallelic\n (n=94)") + theme(legend.position = "none") 
b <- PrintGraph(MTR_AUC, "DD Monoallelic\n (n=435)") + theme(legend.position = "none") + labs(y = "")
c <- PrintGraph(MTR_AUC, "ASD Monoallelic\n (n=190)") + theme(legend.position = "none") + labs(y = "")
d <- PrintGraph(MTR_AUC, "Haploinsufficient\n (n=360)") + theme(legend.position = "none") + labs(y = "")
e <- PrintGraph(MTR_AUC, "Mouse Essential\n (n=2454)") + theme(legend.position = "none") + labs(y = "")

library(patchwork)

patch <- (a|b|c|d|e)

patch

# Saving Figure as PDF
ggsave("Figure2_MTR.pdf", plot = patch, width = 174, height = 87, units = "mm")
