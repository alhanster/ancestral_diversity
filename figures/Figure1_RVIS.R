################################################################################
# Purpose: Figure 1: gnomAD and UKB RVIS                                       #
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

# Pull Functions for gnomAD RVIS Computation, Logistic Regression, DeLong Test
source("scripts/gnomAD_RVIS.R")

# Loading gnomAD Data
gnomAD_data <- fread("data/gnomAD_RVIS.csv")

# Figure 1A: gnomAD RVIS Curves
df.tallies <- gnomAD_data |>
  select(-x, -rvis_afr, -rvis_amr, -rvis_eas,
         -rvis_asj, -rvis_fin, -rvis_nfe, -rvis_sas)


# Calculate RVIS for each ancestry
gnomAD_RVIS <- df.tallies |> 
    mutate( 
      rvis_afr = CalcRVIS(df.tallies, "afr_y"),
      rvis_eas = CalcRVIS(df.tallies, "eas_y"),
      rvis_asj = CalcRVIS(df.tallies, "asj_y"),
      rvis_nfe = CalcRVIS(df.tallies, "nfe_y"),
      rvis_sas = CalcRVIS(df.tallies, "sas_y"),
      rvis_amr = CalcRVIS(df.tallies, "amr_y"),
      rvis_fin = CalcRVIS(df.tallies, "fin_y")
    )

write.csv(gnomAD_RVIS, "output/gnomAD_RVIS_Score.csv", row.names = FALSE)

xy <- df.tallies |> 
  pivot_longer(cols = -c(Gene, mutability), names_to = "ancestry", values_to = "y") %>% 
  filter(ancestry %in% c("afr_y", "sas_y", "amr_y","eas_y", "asj_y", "nfe_y", "fin_y")) %>% 
  mutate(ancestry = toupper(gsub("_y", "", ancestry))) %>% 
  rename(Ancestry = "ancestry")

xy$Ancestry <- factor(xy$Ancestry, levels = c("AFR", "SAS", "AMR", "EAS", "ASJ", "NFE", "FIN"))

figure_1a <- xy |> 
  ggplot(aes(x=mutability, y=y, col=Ancestry)) +
  geom_point(size = 0.5, alpha = 0.1) + 
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  ylim(c(0, 400)) +
  xlim(c(0, 0.0005)) +
  coord_cartesian(xlim = c(0, 0.0003), ylim = c(0, 120)) + 
  xlab("Mutability") + 
  ylab("Common (MAF>0.05%) \n functional variants")+
  scale_color_manual(values = colorblind_palette)+
  theme_classic(base_size=7)+
  theme(
    text = element_text(size = 7),  # Ensure all text is Arial and size 7
    axis.title = element_text(size = 7),  # Axis titles
    axis.text = element_text(size = 7),  # Axis labels
    legend.title = element_text(size = 7),  # Legend title
    legend.text = element_text(size = 7),  # Legend labels
    legend.key.size = unit(3, "mm")  # Legend key size
  )

# Logistic Regression and DeLong Test
dee_results <- PrintAUC(gnomAD_data, score_cols, dee_monoallelic)
dee_AUC <- dee_results$auc_df %>% 
  rename(`DEE Monoallelic\n (n=94)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dee_delong <- dee_results$delong_results %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

dd_results <- PrintAUC(gnomAD_data, score_cols, dd_monoallelic)
dd_AUC <- dd_results$auc_df %>% 
  rename(`DD Monoallelic\n (n=435)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dd_delong <- dd_results$delong_results %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

asd_results <- PrintAUC(gnomAD_data, score_cols, asd_monoallelic)
asd_AUC <- asd_results$auc_df %>% 
  rename(`ASD Monoallelic\n (n=190)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
asd_delong <- asd_results$delong_results %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

mgi_results <- PrintAUC(gnomAD_data, score_cols, mgi_essential)
mgi_AUC <- mgi_results$auc_df %>% 
  rename(`Mouse Essential\n (n=2454)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
mgi_delong <- mgi_results$delong_results %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

HI_results <- PrintAUC(gnomAD_data, score_cols, clingen_HI)
HI_AUC <- HI_results$auc_df %>% 
  rename(`Haploinsufficient\n (n=390)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
HI_delong <- HI_results$delong_results %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

# Compiling AUC Values for Logistic Regression
gnomAD_RVIS_AUC <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
gnomAD_RVIS_AUC <- gnomAD_RVIS_AUC %>% 
  rename("Ancestry" = score) %>% 
  mutate(Ancestry = toupper(gsub("rvis_", "", Ancestry)))

gnomAD_RVIS_AUC$Ancestry <- factor(gnomAD_RVIS_AUC$Ancestry, levels = c("AFR", "SAS", "AMR", "EAS", "ASJ", "NFE", "FIN"))
gnomAD_RVIS_AUC$`Gene List` <- factor(gnomAD_RVIS_AUC$`Gene List`, 
                                      levels = c("DEE Monoallelic\n (n=94)", "DD Monoallelic\n (n=435)", "ASD Monoallelic\n (n=190)", "Haploinsufficient\n (n=390)", "Mouse Essential\n (n=2454)"))


# Figure 1B: gnomAD RVIS Performances by Ancestry
figure_1b <- ggplot(gnomAD_RVIS_AUC, aes(x = `Gene List`, y = AUC, color = Ancestry)) +
  geom_point(position = position_jitter(width = 0.25), size = 1, alpha = 0.7) +
  labs(y = "AUC Scores") +
  ylim(0.6, 0.95)+
  theme_classic(base_size=7)+
  theme(
    text = element_text(size = 7),  # Ensure all text is Arial and size 7
    axis.title = element_text(size = 7),  # Axis titles
    axis.text = element_text(size = 7),  # Axis labels
    legend.title = element_text(size = 7),  # Legend title
    legend.text = element_text(size = 7),  # Legend labels
    legend.key.size = unit(3, "mm"),  # Legend key size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Adjust x-axis text angle & size
    axis.title.x = element_blank()  # Remove x-axis title
  ) +
  scale_color_manual(values = colorblind_palette)


# Compiling Delong Test
delongtest <- rbind(dee_delong, dd_delong, asd_delong, mgi_delong, HI_delong)
delongtest$Adjusted.P.val <- p.adjust(delongtest$De.Long.P.val, method = "BH")

write.csv(delongtest, "output/gnomAD_RVIS_DeLongTest.csv", row.names = FALSE)

# Logistic Regression
dee_AUC <- PrintLogRegResults(gnomAD_data, score_cols, dee_monoallelic) %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

dd_AUC <- PrintLogRegResults(gnomAD_data, score_cols, dd_monoallelic) %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

asd_AUC <- PrintLogRegResults(gnomAD_data, score_cols, asd_monoallelic) %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

mgi_AUC <- PrintLogRegResults(gnomAD_data, score_cols, mgi_essential) %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

HI_AUC <- PrintLogRegResults(gnomAD_data, score_cols, clingen_HI) %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

# Compiling Logistic Regression
gnomAD_RVIS_log <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
gnomAD_RVIS_log <- gnomAD_RVIS_log %>% 
  rename("Ancestry" = score)

write.csv(gnomAD_RVIS_log, "output/gnomAD_RVIS_LogRegression.csv", row.names = FALSE)

# Pull Functions for UKB RVIS Computation, Logistic Regression, DeLong Test
UKB_data <- fread("data/UKB_RVIS.csv")

# Loading UKBiobank Data
source("scripts/UKB_RVIS.R")


# Figure 1C: UKBiobank RVIS Curves
df.tallies <- UKB_data |>
  select(-x, -rvis_afr, -rvis_eas,
         -rvis_asj, -rvis_nfe, -rvis_sas)

UKB_RVIS <- df.tallies |> 
    mutate( 
      rvis_afr = CalcRVIS(df.tallies, "afr_y"),
      rvis_eas = CalcRVIS(df.tallies, "eas_y"),
      rvis_asj = CalcRVIS(df.tallies, "asj_y"),
      rvis_nfe = CalcRVIS(df.tallies, "nfe_y"),
      rvis_sas = CalcRVIS(df.tallies, "sas_y"))

write.csv(UKB_RVIS, "output/UKB_RVIS_Score.csv", row.names = FALSE)

xy <- df.tallies |> 
  pivot_longer(cols = -c(Gene, mutability), names_to = "ancestry", values_to = "y")%>% 
  mutate(ancestry = toupper(gsub("_y", "", ancestry))) %>% 
  rename(Ancestry = "ancestry")

xy$Ancestry <- factor(xy$Ancestry, levels = c("AFR", "SAS", "EAS", "ASJ", "NFE"))

figure_1c <- xy |> 
  ggplot(aes(x=mutability, y=y, col=Ancestry)) +
  geom_point(size = 0.5, alpha = 0.1) + 
  geom_smooth(method = "lm", se = FALSE, size = 0.5) + 
  ylim(c(0, 400)) +
  xlim(c(0, 0.0005)) +
  coord_cartesian(xlim = c(0, 0.0003), ylim = c(0, 120)) + 
  xlab("Mutability") + 
  ylab("Common (MAF>0.05%) \n functional variants")+
  scale_color_manual(values = ukb_colorblind_palette)+
  theme_classic(base_size=7)+
  theme(
    text = element_text(size = 7),  # Ensure all text is Arial and size 7
    axis.title = element_text(size = 7),  # Axis titles
    axis.text = element_text(size = 7),  # Axis labels
    legend.title = element_text(size = 7),  # Legend title
    legend.text = element_text(size = 7),  # Legend labels
    legend.key.size = unit(3, "mm")  # Legend key size
  )

# Logistic Regression and DeLong Test
dee_results <- PrintAUC(UKB_data, score_cols, dee_monoallelic)
dee_AUC <- dee_results$auc_df %>% 
  rename(`DEE Monoallelic\n (n=94)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dee_delong <- dee_results$delong_results %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

dd_results <- PrintAUC(UKB_data, score_cols, dd_monoallelic)
dd_AUC <- dd_results$auc_df %>% 
  rename(`DD Monoallelic\n (n=435)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
dd_delong <- dd_results$delong_results %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

asd_results <- PrintAUC(UKB_data, score_cols, asd_monoallelic)
asd_AUC <- asd_results$auc_df %>% 
  rename(`ASD Monoallelic\n (n=190)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
asd_delong <- asd_results$delong_results %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

mgi_results <- PrintAUC(UKB_data, score_cols, mgi_essential)
mgi_AUC <- mgi_results$auc_df %>% 
  rename(`Mouse Essential\n (n=2454)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
mgi_delong <- mgi_results$delong_results %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

HI_results <- PrintAUC(UKB_data, score_cols, clingen_HI)
HI_AUC <- HI_results$auc_df %>% 
  rename(`Haploinsufficient\n (n=390)` = AUC) %>% 
  pivot_longer(cols = -1, names_to = "Gene List", values_to = "AUC")
HI_delong <- HI_results$delong_results %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(`Gene List`, Score1, Score2, De.Long.P.val)

# Compiling AUC Values for Logistic Regression
UKB_RVIS_AUC <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
UKB_RVIS_AUC <- UKB_RVIS_AUC %>% 
  rename("Ancestry" = score) %>% 
  mutate(Ancestry = toupper(gsub("rvis_", "", Ancestry)))

UKB_RVIS_AUC$Ancestry <- factor(UKB_RVIS_AUC$Ancestry, levels = c("AFR", "SAS", "EAS", "ASJ", "NFE"))
UKB_RVIS_AUC$`Gene List` <- factor(UKB_RVIS_AUC$`Gene List`, 
                                         levels = c("DEE Monoallelic\n (n=94)", "DD Monoallelic\n (n=435)", "ASD Monoallelic\n (n=190)", "Haploinsufficient\n (n=390)", "Mouse Essential\n (n=2454)"))


# Figure 1D: UKBiobank RVIS Performances by Ancestry
figure_1d <- ggplot(UKB_RVIS_AUC, aes(x = `Gene List`, y = AUC, color = Ancestry)) +
  geom_point(position = position_jitter(width = 0.25), size = 1, alpha = 0.7) +
  labs(y = "AUC Scores") +
  scale_fill_manual(values = c("Ancestry A" = "blue", "Ancestry B" = "red")) +
  ylim(0.6, 0.95)+
  theme_classic(base_size=7)+
  theme(
    text = element_text(size = 7),  # Ensure all text is Arial and size 7
    axis.title = element_text(size = 7),  # Axis titles
    axis.text = element_text(size = 7),  # Axis labels
    legend.title = element_text(size = 7),  # Legend title
    legend.text = element_text(size = 7),  # Legend labels
    legend.key.size = unit(3, "mm"),  # Legend key size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Adjust x-axis text angle & size
    axis.title.x = element_blank()  # Remove x-axis title
  )+
  scale_color_manual(values = ukb_colorblind_palette)

# Compiling Delong Test
delongtest <- rbind(dee_delong, dd_delong, asd_delong, mgi_delong, HI_delong)
delongtest$Adjusted.P.val <- p.adjust(delongtest$De.Long.P.val, method = "BH")

write.csv(delongtest, "output/UKB_RVIS_DeLongTest.csv", row.names = FALSE)

# Logistic Regression
dee_AUC <- PrintLogRegResults(UKB_data, score_cols, dee_monoallelic) %>% 
  mutate(`Gene List` = "DEE Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

dd_AUC <- PrintLogRegResults(UKB_data, score_cols, dd_monoallelic) %>% 
  mutate(`Gene List` = "DD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

asd_AUC <- PrintLogRegResults(UKB_data, score_cols, asd_monoallelic) %>% 
  mutate(`Gene List` = "ASD Monoallelic") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

mgi_AUC <- PrintLogRegResults(UKB_data, score_cols, mgi_essential) %>% 
  mutate(`Gene List` = "Mouse Essential") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

HI_AUC <- PrintLogRegResults(UKB_data, score_cols, clingen_HI) %>% 
  mutate(`Gene List` = "Haploinsufficient") %>% 
  select(score, `Gene List`, Log.Reg.P.val, AUC)

# Compiling Logisitic Regression
UKB_RVIS_log <- rbind(dee_AUC, dd_AUC, asd_AUC, mgi_AUC, HI_AUC)
UKB_RVIS_log <- UKB_RVIS_log %>% 
  rename("Ancestry" = score)

write.csv(UKB_RVIS_log, "output/UKB_RVIS_LogRegression.csv", row.names = FALSE)

# Compiling Figures Together
library(patchwork)

patch <- (figure_1a + theme(axis.title.x = element_text(margin = margin(t = -10, unit = "mm")))| figure_1b) / (figure_1c + theme(axis.title.x = element_text(margin = margin(t = -10, unit = "mm"))) | figure_1d) + 
  plot_annotation(tag_levels = 'a', theme = theme(plot.tag = element_text(size = 7)))

# Save Figure as PDF
ggsave("figure1.pdf", plot = patch, path = "output", width = 180, height = 160, units = "mm")
