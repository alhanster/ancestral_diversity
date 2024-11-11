################################################################################
# Purpose: Supplementary Figure 3: UKB LOF FDR                                 #
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

# Loading UKB Data for LOF FDR
data <- fread("data/UKB_LOF_FDR.csv")

df <- data |>
  select(-LOF_FDR_AFR, -LOF_FDR_ASJ, -LOF_FDR_EAS, -LOF_FDR_SAS, -LOF_FDR_Diverse43k, 
         -LOF_FDR_NFE20k, -LOF_FDR_NFE43k, -LOF_FDR_NFE440k, -LOF_FDR_FullDataset)

# Computing LOF FDR
afr_lof_fdr <- df %>% 
  filter(afr_total != 0) %>% 
  rowwise() %>%
  mutate(afr_pval = binom.test(x = afr_lof, n = afr_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, afr_lof, afr_total, afr_pval)

sas_lof_fdr <- df %>% 
  filter(sas_total != 0) %>% 
  rowwise() %>%
  mutate(sas_pval = binom.test(x = sas_lof, n = sas_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, sas_lof, sas_total, sas_pval)

asj_lof_fdr <- df %>% 
  filter(asj_total != 0) %>% 
  rowwise() %>%
  mutate(asj_pval = binom.test(x = asj_lof, n = asj_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, asj_lof, asj_total, asj_pval)

eas_lof_fdr <- df %>% 
  filter(eas_total != 0) %>% 
  rowwise() %>%
  mutate(eas_pval = binom.test(x = eas_lof, n = eas_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, eas_lof, eas_total, eas_pval)

max_diverse_lof_fdr <- df %>% 
  filter(max_diverse_total != 0) %>% 
  rowwise() %>%
  mutate(max_diverse_pval = binom.test(x = max_diverse_lof, n = max_diverse_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, max_diverse_lof, max_diverse_total, max_diverse_pval)

nfe_20k_lof_fdr <- df %>% 
  filter(nfe_20k_total != 0) %>% 
  rowwise() %>%
  mutate(nfe_20k_pval = binom.test(x = nfe_20k_lof, n = nfe_20k_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, nfe_20k_lof, nfe_20k_total, nfe_20k_pval)

nfe_43k_lof_fdr <- df %>% 
  filter(nfe_43k_total != 0) %>% 
  rowwise() %>%
  mutate(nfe_43k_pval = binom.test(x = nfe_43k_lof, n = nfe_43k_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, nfe_43k_lof, nfe_43k_total, nfe_43k_pval)

nfe_440k_lof_fdr <- df %>% 
  filter(nfe_440k_total != 0) %>% 
  rowwise() %>%
  mutate(nfe_440k_pval = binom.test(x = nfe_440k_lof, n = nfe_440k_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, nfe_440k_lof, nfe_440k_total, nfe_440k_pval)

all_lof_fdr <- df %>% 
  filter(all_total != 0) %>% 
  rowwise() %>%
  mutate(all_pval = binom.test(x = all_lof, n = all_total, p = exp_lof_percent, alternative = "less")$p.value) %>% 
  select(Gene, all_lof, all_total, all_pval)

lof.fdr <- all_lof_fdr %>% 
  full_join(afr_lof_fdr, by = "Gene") %>% 
  full_join(sas_lof_fdr, by = "Gene") %>% 
  full_join(asj_lof_fdr, by = "Gene") %>% 
  full_join(eas_lof_fdr, by = "Gene") %>% 
  full_join(max_diverse_lof_fdr, by = "Gene") %>% 
  full_join(nfe_20k_lof_fdr, by = "Gene") %>% 
  full_join(nfe_43k_lof_fdr, by = "Gene") %>% 
  full_join(nfe_440k_lof_fdr, by = "Gene") %>% 
  rename(
    AFR = afr_pval,
    ASJ = asj_pval,
    EAS = eas_pval,
    SAS = sas_pval,
    `Maximally Diverse (n=43k)` = max_diverse_pval,
    `NFE (n=20k)` = nfe_20k_pval,
    `NFE (n=43k)` = nfe_43k_pval,
    `NFE (n=440k)` = nfe_440k_pval,
    `Full Dataset (n=460k)` = all_pval
  ) %>% 
  select(Gene, `Maximally Diverse (n=43k)`, `NFE (n=43k)`, `NFE (n=440k)`, `Full Dataset (n=460k)`, AFR, ASJ, EAS, SAS, `NFE (n=20k)`)

write.csv(lof.fdr, "output/LOF_FDR_Score.csv", row.names = FALSE)

# Logistic Regression and DeLong Test
dee_results <- PrintAUC(lof.fdr, score_cols, dee_monoallelic)
dee_AUC <- dee_results$auc_df %>% 
  rename(`DEE Monoallelic\n (n=94)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dee_delong <- dee_results$delong_results %>% 
  mutate(Gene_List = "DEE Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

dd_results <- PrintAUC(lof.fdr, score_cols, dd_monoallelic)
dd_AUC <- dd_results$auc_df %>% 
  rename(`DD Monoallelic\n (n=435)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dd_delong <- dd_results$delong_results %>% 
  mutate(Gene_List = "DD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

asd_results <- PrintAUC(lof.fdr, score_cols, asd_monoallelic)
asd_AUC <- asd_results$auc_df %>% 
  rename(`ASD Monoallelic\n (n=190)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
asd_delong <- asd_results$delong_results %>% 
  mutate(Gene_List = "ASD Monoallelic") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

mgi_results <- PrintAUC(lof.fdr, score_cols, mgi_essential)
mgi_AUC <- mgi_results$auc_df %>% 
  rename(`Mouse Essential\n (n=2454)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
mgi_delong <- mgi_results$delong_results %>% 
  mutate(Gene_List = "Mouse Essential") %>% 
  select(Gene_List, Score1, Score2, De.Long.P.val)

HI_results <- PrintAUC(lof.fdr, score_cols, clingen_HI)
HI_AUC <- HI_results$auc_df %>% 
  rename(`Haploinsufficient\n (n=390)` = AUC) %>% 
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
write.csv(delongtest, "output/LOF_FDR_DeLongTest.csv", row.names = FALSE)

# Logistic Regression
dee_AUC <- PrintLogRegResults(lof.fdr, score_cols, dee_monoallelic) %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

dd_AUC <- PrintLogRegResults(lof.fdr, score_cols, dd_monoallelic) %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

asd_AUC <- PrintLogRegResults(lof.fdr, score_cols, asd_monoallelic) %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

mgi_AUC <- PrintLogRegResults(lof.fdr, score_cols, mgi_essential) %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

HI_AUC <- PrintLogRegResults(lof.fdr, score_cols, clingen_HI) %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)


# Compiling Logistic Regression
UKB_LOF_FDR_log <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
UKB_LOF_FDR_log <- UKB_LOF_FDR_log %>% 
  rename("Ancestry" = score)
write.csv(UKB_LOF_FDR_log, "output/LOF_FDR_LogRegression.csv", row.names = FALSE)


# Figure 4B: UKB LOF FDR Performances by Ancestry
a <- PrintGraph(LOF_FDR_AUC, "DEE Monoallelic\n (n=94)") + theme(legend.position = "none") + labs(y = "LOF-FDR AUC Scores")
b <- PrintGraph(LOF_FDR_AUC, "DD Monoallelic\n (n=435)") + theme(legend.position = "none") + labs(y = "")
c <- PrintGraph(LOF_FDR_AUC, "ASD Monoallelic\n (n=190)") + theme(legend.position = "none") + labs(y = "")
d <- PrintGraph(LOF_FDR_AUC, "Haploinsufficient\n (n=390)") + theme(legend.position = "none") + labs(y = "")
e <- PrintGraph(LOF_FDR_AUC, "Mouse Essential\n (n=2454)") + theme(legend.position = "none") + labs(y = "")

library(patchwork)
patch <- (a|b|c|d|e)

# Saving Figures as PDFs
ggsave("Supp_Figure_LOF_FDR.pdf", plot = patch, path = "output", width = 174, height = 87, units = "mm")
