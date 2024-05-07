################################################################################
# Purpose:                                                                     #
# Author: Ryan Dhindsa                                                         #
################################################################################

# Imports ----------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(pROC)
library(cowplot)
library(PRROC)

# Globals ----------------------------------------------------------------------
# TODO: make sure these are complete
SYN <- c("synonymous_variant")
MISSENSE <- c("missense_variant", "inframe_deletion", "inframe_insertion")
LOF <- c("frameshift_variant", "stop_gained", "start_lost", 
       "splice_acceptor_variant", "splice_donor_variant",
       "stop_gained&frameshift_variant", "splice_donor_variant&coding_sequence_variant&intron_variant", 
       "splice_donor_variant&intron_variant", "stop_gained&inframe_insertion", "frameshift_variant&stop_lost", 
       "frameshift_variant&start_lost", "start_lost&splice_region_variant", 
       "stop_gained&protein_altering_variant", "stop_gained&frameshift_variant&splice_region_variant", 
       "stop_gained&inframe_deletion", "frameshift_variant&stop_retained_variant")

# Functions --------------------------------------------------------------------
CalcXY <- function(df, maf.cutoff) {
  df.x <- df |> 
    group_by(Gene) |> 
    summarise(x = n()) 
  
  df.y <- df |> 
    filter(Consequence %in% c(MISSENSE, LOF)) |> 
    group_by(Gene) |>
    summarise(
      afr_y = sum(AF_afr > maf.cutoff, na.rm = T),
      ami_y = sum(AF_ami > maf.cutoff, na.rm = T),
      amr_y = sum(AF_amr > maf.cutoff, na.rm = T), 
      eas_y = sum(AF_eas > maf.cutoff, na.rm = T), 
      asj_y = sum(AF_asj > maf.cutoff, na.rm = T), 
      fin_y = sum(AF_fin > maf.cutoff, na.rm = T), 
      nfe_y = sum(AF_nfe > maf.cutoff, na.rm = T), 
      mid_y = sum(AF_mid > maf.cutoff, na.rm = T), 
      sas_y = sum(AF_sas > maf.cutoff, na.rm = T),
      oth_y = sum(AF_oth > maf.cutoff, na.rm = T), 
      all_y = sum(AF > maf.cutoff, na.rm = T),
      popmax_y = sum(AF_popmax > maf.cutoff, na.rm = T)
      )
  
  tallies <- df.x |> 
    left_join(df.y)
  
  return(tallies)
}


CalcRVIS <- function(df, y.colname) {
  df <- df |> 
    replace(is.na(df), 0)
  
  fm <- as.formula(paste(y.colname, "~", "x"))
  model <- lm(fm, data = df)
  rvis.scores <- rstudent(model)
  
  return(rvis.scores)
}



PlotROCs <- function(df, score_cols, gene.list, list.name="Test") {
  df_long <- df %>% 
    select(c(Gene, tidyselect::all_of(score_cols))) |> 
    pivot_longer(cols = score_cols, 
                 names_to = "score_type", 
                 values_to = "score_value")
  
  df_long <- df_long |> 
    mutate(gene_list = ifelse(Gene %in% gene.list$V1, 1, 0))
  
  # Compute the ROC curve and AUC for each score
  roc_list <- lapply(unique(df_long$score_type), function(col) {
    df_sub <- df_long %>% filter(score_type == col)
    
    roc_obj <- roc(df_sub$gene_list ~ df_sub$score_value)
    auc_val <- auc(roc_obj)
    
    # Create data frame for AUC
    auc_df <- data.frame(score = col, AUC = auc_val)
    
    # Create data frame for ROC curve
    roc_df <- data.frame(score = col, FPR = 1 - roc_obj$sensitivities, TPR = roc_obj$specificities)
    
    list(auc_df = auc_df, roc_df = roc_df)
  })
  
  # Combine all ROC curves and AUCs
  auc_df <- do.call(rbind, lapply(roc_list, `[[`, "auc_df"))
  roc_df <- do.call(rbind, lapply(roc_list, `[[`, "roc_df"))
  
  # Print AUCs
  print(auc_df)
  
  # Plot the ROC curves
  roc_df <- roc_df |> 
    left_join(auc_df) |> 
    mutate(score = paste0(score, " (AUC=", signif(AUC, 2), ")")) |> 
    arrange(AUC)
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR, color = reorder(score, AUC, decreasing = T))) +
    geom_line() +
    theme_minimal() +
    ggtitle(list.name) +
    labs(x = "False Positive Rate", y = "True Positive Rate", color = "Score type") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    ggthemes::scale_color_tableau()
  
  p
}


PlotPRCurves <- function(df, score_cols, gene.list, list.name = "Test") {
  df_long <- df %>%
    select(c(Gene, all_of(score_cols))) |>
    pivot_longer(cols = score_cols,
                 names_to = "score_type",
                 values_to = "score_value")
  
  df_long <- df_long |>
    mutate(gene_list = ifelse(Gene %in% gene.list$V1, 1, 0))
  
  # Compute the precision-recall curve and AUC for each score
  pr_list <- lapply(unique(df_long$score_type), function(col) {
    df_sub <- df_long %>% filter(score_type == col)
    
    pr_obj <- roc(df_sub$gene_list, df_sub$score_value, direction = "<")
    auc_val <- auc(pr_obj, method = "PR")
    
    # Create data frame for AUC
    auc_df <- data.frame(score = col, AUC = auc_val)
    
    # Create data frame for precision-recall curve
    pr_df <- data.frame(score = col, Precision = pr_obj$specificities, Recall = pr_obj$sensitivities)
    
    list(auc_df = auc_df, pr_df = pr_df)
  })
  
  # Combine all precision-recall curves and AUCs
  auc_df <- do.call(rbind, lapply(pr_list, `[[`, "auc_df"))
  pr_df <- do.call(rbind, lapply(pr_list, `[[`, "pr_df"))
  
  # Print AUCs
  print(auc_df)
  
  # Plot the precision-recall curves
  pr_df <- pr_df |>
    left_join(auc_df) |>
    mutate(score = paste0(score, " (AUC=", signif(AUC, 2), ")")) |>
    arrange(AUC)
  
  p <- ggplot(pr_df, aes(x = Recall, y = Precision, color = reorder(score, AUC, decreasing = TRUE))) +
    geom_line() +
    theme_minimal() +
    ggtitle(list.name) +
    labs(x = "Recall", y = "Precision", color = "Score type") +
    geom_segment(aes(x = 0, y = max(pr_df$Precision), xend = 1, yend = max(pr_df$Precision)), linetype = "dashed", color = "grey") +
    geom_segment(aes(x = max(pr_df$Recall), y = 0, xend = max(pr_df$Recall), yend = 1), linetype = "dashed", color = "grey") +
    geom_point(data = pr_df %>% filter(Recall == max(Recall) | Precision == max(Precision)), size = 3, shape = 21, fill = "white") +
    geom_text(data = pr_df %>% filter(Recall == max(Recall) | Precision == max(Precision)), aes(label = score), nudge_x = 0.02, nudge_y = 0.02, size = 3) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggthemes::scale_color_tableau()
  
  return(p)
}



# # Main -------------------------------------------------------------------------
# setwd("~/Desktop/Dhindsa Lab/RVIS/")
# gnomad.metrics <- fread("Data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
# # gnomad.maf <- fread("Data/gnomad_exomes_2.1.1.tsv.gz")
# 
# gnomad.filtered <- gnomad.maf |> 
#   filter(Transcript %in% gnomad.metrics$transcript) |> 
#   filter(Consequence %in% c(SYN, MISSENSE, LOF)) |> 
#   filter(FILTER %in% c("None"))
# 
# 
# # Calculate X and Y for RVIS calculation
# df.tallies <- gnomad.filtered |> 
#   CalcXY(maf.cutoff = 0.001)
# 
# # Calculate RVIS for each ancestry
# df <- df.tallies |> 
#   mutate(rvis_all = CalcRVIS(df.tallies, "all_y"), 
#          rvis_popmax = CalcRVIS(df.tallies, "popmax_y"), 
#          rvis_afr = CalcRVIS(df.tallies, "afr_y"),
#          rvis_amr = CalcRVIS(df.tallies, "amr_y"),
#          rvis_eas = CalcRVIS(df.tallies, "eas_y"),
#          rvis_asj = CalcRVIS(df.tallies, "asj_y"),
#          rvis_fin = CalcRVIS(df.tallies, "fin_y"),
#          rvis_nfe = CalcRVIS(df.tallies, "nfe_y"),
#          rvis_sas = CalcRVIS(df.tallies, "sas_y"),
#          rvis_ami = CalcRVIS(df.tallies, "ami_y"),
#          rvis_mid = CalcRVIS(df.tallies, "mid_y"),
#          rvis_oth = CalcRVIS(df.tallies, "oth_y"))


# cor.res <- cor(df[,11:18])

# Plot XY tallies

# xy <- df.tallies |> 
#   pivot_longer(cols = -c(Gene, x), names_to = "ancestry", values_to = "y") 
# 
# xy |> 
#   filter(ancestry != "popmax_y") |> 
#   ggplot(aes(x=x, y=y, col=ancestry)) +
#   geom_point(size = 0.5, alpha = 0.05) + 
#   geom_smooth(method = "lm") + 
#   coord_cartesian(xlim = c(0, 5000), ylim = c(0, 150)) + 
#   ggthemes::scale_color_tableau() +
#   theme_bw() + 
#   xlab("All observed variants") + 
#   ylab("Common (MAF>0.1%) functional variants")
