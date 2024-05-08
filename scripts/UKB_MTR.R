################################################################################
# Purpose: UKB MTR                                                             #
# Author: Alexander Han                                                        #
################################################################################


score_cols = c("Maximally Diverse (n=43k)", "NFE (n=43k)", "NFE (n=440k)", 
               "Full Dataset (n=460k)", "AFR", "ASJ", "EAS", "SAS", "NFE (n=20k)")


PrintAUC <- function(MTR, score_cols, gene.list) {
  
  df_long <- MTR %>% 
    select(c(Gene, all_of(score_cols))) |> 
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
    
    list(roc_obj = roc_obj, auc_df = auc_df, roc_df = roc_df)
  })
  
  # Combine all ROC curves and AUCs
  auc_df <- do.call(rbind, lapply(roc_list, `[[`, "auc_df"))
  roc_df <- do.call(rbind, lapply(roc_list, `[[`, "roc_df"))
  
  # Sort AUCs in descending order
  auc_df <- auc_df %>% 
    arrange(desc(AUC))
  
  # Perform DeLong test for each pair of ROC curves and store results in a new table
  delong_results <- data.frame()
  for (i in 1:(length(roc_list) - 1)) {
    for (j in (i + 1):length(roc_list)) {
      test_result <- roc.test(roc_list[[i]]$roc_obj, roc_list[[j]]$roc_obj)
      delong_results <- rbind(delong_results, data.frame(
        Score1 = roc_list[[i]]$auc_df$score,
        Score2 = roc_list[[j]]$auc_df$score,
        De.Long.P.val = test_result$p.value
      ))
    }
  }
  
  return(list(auc_df = auc_df, roc_df = roc_df, delong_results = delong_results))
}


PrintLogRegResults <- function(MTR, score_cols, gene.list, maf.threshold = 0.0005) {
  
  df_long <- MTR %>% 
    select(c(Gene, all_of(score_cols))) |> 
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
    
    list(roc_obj = roc_obj, auc_df = auc_df, roc_df = roc_df)
  })
  
  # Combine all ROC curves and AUCs
  auc_df <- do.call(rbind, lapply(roc_list, `[[`, "auc_df"))
  roc_df <- do.call(rbind, lapply(roc_list, `[[`, "roc_df"))
  
  # Sort AUCs in descending order
  auc_df <- auc_df %>% 
    arrange(desc(AUC))
  
  # Logistic Regression Analysis
  log_reg_results <- lapply(unique(df_long$score_type), function(col) {
    df_sub <- df_long %>% filter(score_type == col)
    
    model <- glm(gene_list ~ score_value, data = df_sub, family = binomial)
    summary_df <- summary(model)$coefficients
    
    # Create a dataframe to store results
    results_df <- data.frame(score = col, Estimate = summary_df["score_value", "Estimate"],
                             Std.Error = summary_df["score_value", "Std. Error"],
                             z.value = summary_df["score_value", "z value"],
                             Log.Reg.P.val = summary_df["score_value", "Pr(>|z|)"])
    
    return(results_df)
  })
  
  log_reg_df <- do.call(rbind, log_reg_results)
  
  log_reg_df <- log_reg_df %>%
    select(score, Log.Reg.P.val) %>% 
    left_join(auc_df, by = c("score"))
  
  # Print Logistic Regression Results
  return(log_reg_df)
}

colorblind_palette <- c(
  "#56B4E9", # sky blue
  "#E69F00", # orange
  "#009E73", # bluish green
  "#CC79A7", # reddish purple
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#999999", # light grey
  "#66C2A5"  # turquoise
)

PrintGraph <- function(MTR_AUC, gene.list) {
  
  MTR_AUC_gene <- MTR_AUC %>% 
    filter(`Gene List` %in% paste(gene.list)) %>% 
    arrange(AUC)
  
  MTR_AUC_ancestry <- MTR_AUC_gene %>% 
    pull(Ancestry)
  
  MTR_AUC_gene$Ancestry <- factor(MTR_AUC_gene$Ancestry, levels = MTR_AUC_ancestry)
  
  figure <- ggplot(MTR_AUC_gene, aes(x = Ancestry, y = AUC, color = Ancestry)) +
    geom_point(shape = 15, size = 2) +
    labs(y = "Gene-level MTR AUC Scores") +
    scale_color_manual(values = c(
      "ASJ" = "#56B4E9", 
      "EAS" = "#E69F00",
      "SAS" = "#009E73", 
      "AFR" = "#CC79A7", 
      "NFE (n=20k)" = "#F0E442", 
      "NFE (n=43k)" = "#0072B2", 
      "NFE (n=440k)"= "#D55E00", 
      "Maximally Diverse (n=43k)" = "#999999", 
      "Full Dataset (n=460k)" = "#66C2A5")
    )+
    ggtitle(paste(gene.list))+
    ylim(0.5,0.9)+
    theme_classic(base_size=7)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5))
  
  figure
}
