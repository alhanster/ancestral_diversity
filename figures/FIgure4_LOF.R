################################################################################
# Purpose: Figure 4: UKB LOF O/E and FDR                                       #
# Author: Alexander Han                                                        #
################################################################################

# Paste Path to Repository
setwd()

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
    `Maximally Diverse (n=43k)` = Diverse_LOF / possible_lof,
    `NFE (n=43k)` = Nfe_43k_LOF / possible_lof,
    `NFE (n=440k)` = Nfe_440k_LOF / possible_lof,
    `Full Dataset (n=460k)` = All_LOF / possible_lof,
    AFR = afr_LOF / possible_lof,
    ASJ = asj_LOF / possible_lof,
    EAS = eas_LOF / possible_lof,
    SAS = sas_LOF / possible_lof,
    `NFE (n=20k)` = nfe_20k_LOF / possible_lof
  ) %>%
  select(Gene, `Maximally Diverse (n=43k)`, `NFE (n=43k)`, `NFE (n=440k)`, `Full Dataset (n=460k)`, 
         AFR, ASJ, EAS, SAS, `NFE (n=20k)`)


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
  rename(`Haploinsufficient\n (n=360)` = AUC) %>% 
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
delongtest

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


# Figure 4a: UKB LOF O/E Performances by Ancestry
a <- PrintGraph(LOF_OE_AUC, "DEE Monoallelic\n (n=94)") + theme(legend.position = "none") 
b <- PrintGraph(LOF_OE_AUC, "DD Monoallelic\n (n=435)") + theme(legend.position = "none") + labs(y = "")
c <- PrintGraph(LOF_OE_AUC, "ASD Monoallelic\n (n=190)") + theme(legend.position = "none") + labs(y = "")
d <- PrintGraph(LOF_OE_AUC, "Haploinsufficient\n (n=360)") + theme(legend.position = "none") + labs(y = "")
e <- PrintGraph(LOF_OE_AUC, "Mouse Essential\n (n=2454)") + theme(legend.position = "none") + labs(y = "")

library(patchwork)

patch1 <- (a|b|c|d|e)


# Loading UKB Data for LOF FDR
data <- fread("data/UKB_LOF_FDR.csv")

df <- data |>
  select(-LOF_FDR_Diverse43k, -LOF_FDR_NFE43k, -LOF_FDR_NFE440k, -LOF_FDR_FullDataset, -LOF_FDR_AFR, -LOF_FDR_ASJ, -LOF_FDR_EAS, -LOF_FDR_SAS, -LOF_FDR_NFE20k)


# Computing LOF FDR
LOF_FDR <- df %>%
  mutate(
    `Maximally Diverse (n=43k)` = (Diverse_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)) ,
    `NFE (n=43k)` = (Nfe_43k_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)),
    `NFE (n=440k)` = (Nfe_440k_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)),
    `Full Dataset (n=460k)` = (All_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)),
    
    AFR = (afr_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)) ,
    ASJ = (asj_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)),
    EAS = (eas_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)),
    SAS = (sas_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis)),
    `NFE (n=20k)` = (nfe_20k_LOF/total_var)/((mu_lof+mu_lof*1.25)/((mu_lof+mu_lof*1.25)+mu_syn+mu_mis))
  ) %>%
  select(Gene, `Maximally Diverse (n=43k)`, `NFE (n=43k)`, `NFE (n=440k)`, `Full Dataset (n=460k)`, AFR, ASJ, EAS, SAS, `NFE (n=20k)`)

# Logistic Regression and DeLong Test
dee_results <- PrintAUC(LOF_FDR, score_cols, dee_monoallelic)
dee_AUC <- dee_results$auc_df %>% 
  rename(`DEE Monoallelic\n (n=94)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dee_delong <- dee_results$delong_results %>% 
  mutate(Gene_List = "DEE Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

dd_results <- PrintAUC(LOF_FDR, score_cols, dd_monoallelic)
dd_AUC <- dd_results$auc_df %>% 
  rename(`DD Monoallelic\n (n=435)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dd_delong <- dd_results$delong_results %>% 
  mutate(Gene_List = "DD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

asd_results <- PrintAUC(LOF_FDR, score_cols, asd_monoallelic)
asd_AUC <- asd_results$auc_df %>% 
  rename(`ASD Monoallelic\n (n=190)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
asd_delong <- asd_results$delong_results %>% 
  mutate(Gene_List = "ASD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

mgi_results <- PrintAUC(LOF_FDR, score_cols, mgi_essential)
mgi_AUC <- mgi_results$auc_df %>% 
  rename(`Mouse Essential\n (n=2454)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
mgi_delong <- mgi_results$delong_results %>% 
  mutate(Gene_List = "Mouse Essential") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

HI_results <- PrintAUC(LOF_FDR, score_cols, clingen_HI)
HI_AUC <- HI_results$auc_df %>% 
  rename(`Haploinsufficient\n (n=360)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
HI_delong <- HI_results$delong_results %>% 
  mutate(Gene_List = "Haploinsufficient") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)


# Compiling AUC Values for Logistic Regression
LOF_FDR_AUC <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
LOF_FDR_AUC <- LOF_FDR_AUC %>% 
  rename("Ancestry" = score)

# Compiling DeLong Test
delongtest <- rbind(dee_delong, dd_delong, asd_delong, mgi_delong, HI_delong)


# Logistic Regression
dee_AUC <- PrintLogRegResults(LOF_FDR, score_cols, dee_monoallelic) %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

dd_AUC <- PrintLogRegResults(LOF_FDR, score_cols, dd_monoallelic) %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

asd_AUC <- PrintLogRegResults(LOF_FDR, score_cols, asd_monoallelic) %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

mgi_AUC <- PrintLogRegResults(LOF_FDR, score_cols, mgi_essential) %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

HI_AUC <- PrintLogRegResults(LOF_FDR, score_cols, clingen_HI) %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)


# Compiling Logistic Regression
UKB_LOF_FDR_log <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
UKB_LOF_FDR_log <- UKB_LOF_FDR_log %>% 
  rename("Ancestry" = score)


# Figure 4B: UKB LOF FDR Performances by Ancestry
a <- PrintGraph(LOF_FDR_AUC, "DEE Monoallelic\n (n=94)") + theme(legend.position = "none") 
b <- PrintGraph(LOF_FDR_AUC, "DD Monoallelic\n (n=435)") + theme(legend.position = "none") + labs(y = "")
c <- PrintGraph(LOF_FDR_AUC, "ASD Monoallelic\n (n=190)") + theme(legend.position = "none") + labs(y = "")
d <- PrintGraph(LOF_FDR_AUC, "Haploinsufficient\n (n=360)") + theme(legend.position = "none") + labs(y = "")
e <- PrintGraph(LOF_FDR_AUC, "Mouse Essential\n (n=2454)") + theme(legend.position = "none") + labs(y = "")

patch2 <- (a|b|c|d|e)

# Compiling Figure
figure4a <- patch1 + plot_annotation('A')
figure4b <- patch2 + plot_annotation('B')


# Saving Figures as PDFs
ggsave("Figure4a.pdf", plot = figure4a, width = 174, height = 87, units = "mm")
ggsave("Figure4b.pdf", plot = figure4b, width = 174, height = 87, units = "mm")
