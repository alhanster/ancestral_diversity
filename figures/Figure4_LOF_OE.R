################################################################################
# Purpose: Figure 4: UKB LOF O/E                                               #
# Author: Alexander Han                                                        #
################################################################################

folder_path <- "output"

# Check if the folder already exists
if (!dir.exists(folder_path)) {
  # Create the folder if it does not exist
  dir.create(folder_path, recursive = TRUE)
}

# Pull Variant Annotation and Gene Lists
source("scripts/GeneList_VariantAnnotations.R")

# Pull Functions for UKB MTR Computation, Logistic Regression, DeLong Test
source("scripts/UKB_LOF.R")

# Loading UKB Data for LOF O/E
data <- fread("data/UKB_LOF_OE.csv")


df <- data |>
  select(-LOF_OE_Diverse43k, -LOF_OE_NFE43k, -LOF_OE_NFE440k, -LOF_OE_FullDataset, -LOF_OE_AFR, -LOF_OE_ASJ, -LOF_OE_EAS, -LOF_OE_SAS, -LOF_OE_NFE20k)


# Computing LOF O/E
LOF_OE <- df %>%
  mutate(
    `Maximally Diverse (n=43k)`= (max_diverse_lof/max_diverse_total)/((mu_lof+mu_lof*1.25)/mutability),
    `NFE (n=43k)` = (nfe_43k_lof/nfe_43k_total)/((mu_lof+mu_lof*1.25)/mutability),
    `NFE (n=440k)` = (nfe_440k_lof/nfe_440k_total)/((mu_lof+mu_lof*1.25)/mutability),
    `Full Dataset (n=460k)` = (all_lof/all_total)/((mu_lof+mu_lof*1.25)/mutability),
    
    `AFR` = (afr_lof/afr_total)/((mu_lof+mu_lof*1.25)/mutability) ,
    `ASJ` = (asj_lof/asj_total)/((mu_lof+mu_lof*1.25)/mutability),
    `EAS` = (eas_lof/eas_total)/((mu_lof+mu_lof*1.25)/mutability),
    `SAS` = (sas_lof/sas_total)/((mu_lof+mu_lof*1.25)/mutability),
    `NFE (n=20k)` = (nfe_20k_lof/nfe_20k_total)/((mu_lof+mu_lof*1.25)/mutability)
  ) %>%
  select(Gene, `Maximally Diverse (n=43k)`, `NFE (n=43k)`, `NFE (n=440k)`, `Full Dataset (n=460k)`, 
         AFR, ASJ, EAS, SAS, `NFE (n=20k)`)

write.csv(LOF_OE, "output/LOF_OE_Score.csv", row.names = FALSE)

# Logistic Regression and DeLong Test
dee_results <- PrintAUC(LOF_OE, score_cols, dee_monoallelic)
dee_AUC <- dee_results$auc_df %>% 
  rename(`DEE Monoallelic\n (n=94)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dee_delong <- dee_results$delong_results %>% 
  mutate(Gene_List = "DEE Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

dd_results <- PrintAUC(LOF_OE, score_cols, dd_monoallelic)
dd_AUC <- dd_results$auc_df %>% 
  rename(`DD Monoallelic\n (n=435)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dd_delong <- dd_results$delong_results %>% 
  mutate(Gene_List = "DD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

asd_results <- PrintAUC(LOF_OE, score_cols, asd_monoallelic)
asd_AUC <- asd_results$auc_df %>% 
  rename(`ASD Monoallelic\n (n=190)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
asd_delong <- asd_results$delong_results %>% 
  mutate(Gene_List = "ASD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

mgi_results <- PrintAUC(LOF_OE, score_cols, mgi_essential)
mgi_AUC <- mgi_results$auc_df %>% 
  rename(`Mouse Essential\n (n=2454)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
mgi_delong <- mgi_results$delong_results %>% 
  mutate(Gene_List = "Mouse Essential") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

HI_results <- PrintAUC(LOF_OE, score_cols, clingen_HI)
HI_AUC <- HI_results$auc_df %>% 
  rename(`Haploinsufficient\n (n=390)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
HI_delong <- HI_results$delong_results %>% 
  mutate(Gene_List = "Haploinsufficient") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

# Compiling AUC Values for Logistic Regression
LOF_OE_AUC <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
LOF_OE_AUC <- LOF_OE_AUC %>% 
  rename("Ancestry" = score)

# Compiling DeLong Test
delongtest <- rbind(dee_delong, dd_delong, asd_delong, mgi_delong, HI_delong)
delongtest$Adjusted.P.val <- p.adjust(delongtest$De.Long.P.val, method = "BH")

write.csv(delongtest, "output/LOF_OE_DeLongTest.csv", row.names = FALSE)

# Logistic Regression
dee_AUC <- PrintLogRegResults(LOF_OE, score_cols, dee_monoallelic) %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

dd_AUC <- PrintLogRegResults(LOF_OE, score_cols, dd_monoallelic) %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

asd_AUC <- PrintLogRegResults(LOF_OE, score_cols, asd_monoallelic) %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

mgi_AUC <- PrintLogRegResults(LOF_OE, score_cols, mgi_essential) %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

HI_AUC <- PrintLogRegResults(LOF_OE, score_cols, clingen_HI) %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

# Compiling Logistic Regression
UKB_LOF_OE_log <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
UKB_LOF_OE_log <- UKB_LOF_OE_log %>% 
  rename("Ancestry" = score)
write.csv(UKB_LOF_OE_log, "output/LOF_OE_LogRegression.csv", row.names = FALSE)


# Figure 4a: UKB LOF O/E Performances by Ancestry
a <- PrintGraph(LOF_OE_AUC, "DEE Monoallelic\n (n=94)") + theme(legend.position = "none") 
b <- PrintGraph(LOF_OE_AUC, "DD Monoallelic\n (n=435)") + theme(legend.position = "none") + labs(y = "")
c <- PrintGraph(LOF_OE_AUC, "ASD Monoallelic\n (n=190)") + theme(legend.position = "none") + labs(y = "")
d <- PrintGraph(LOF_OE_AUC, "Haploinsufficient\n (n=390)") + theme(legend.position = "none") + labs(y = "")
e <- PrintGraph(LOF_OE_AUC, "Mouse Essential\n (n=2454)") + theme(legend.position = "none") + labs(y = "")

library(patchwork)

patch <- (a|b|c|d|e)

ggsave("figure4.pdf", plot = patch, path = "output", width = 180, height = 100, units = "mm")
