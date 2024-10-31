################################################################################
# Purpose: UKB RVIS calculation, Logistic Regression, DeLong Test              #
# Author: Alexander Han                                                        #
################################################################################

score_cols = c("rvis_afr", "rvis_asj", "rvis_eas", "rvis_nfe", "rvis_sas")

ukb_colorblind_palette <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7" # reddish purple
)


CalcRVIS <- function(df, y.colname) {
  
  df <- df |> 
    replace(is.na(df), 0)
  
  fm <- as.formula(paste(y.colname, "~", "mutability"))
  model <- lm(fm, data = df)
  rvis.scores <- rstudent(model)
  
  return(rvis.scores)
}


PrintAUC <- function(df, score_cols, gene.list) {
  
  df.tallies <- df |> 
    select(Gene, x, afr_y, asj_y, eas_y, nfe_y, sas_y, mutability)
  
  # Calculate RVIS for each ancestry
  df <- df.tallies |> 
    mutate( 
      rvis_afr = CalcRVIS(df.tallies, "afr_y"),
      rvis_eas = CalcRVIS(df.tallies, "eas_y"),
      rvis_asj = CalcRVIS(df.tallies, "asj_y"),
      rvis_nfe = CalcRVIS(df.tallies, "nfe_y"),
      rvis_sas = CalcRVIS(df.tallies, "sas_y"))
  
  df_long <- df %>% 
    select(c(Gene, all_of(score_cols))) |> 
    pivot_longer(cols = score_cols, 
                 names_to = "score_type", 
                 values_to = "score_value")|>
    filter(!is.na(score_value))
  
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

PrintLogRegResults <- function(df, score_cols, gene.list) {
  
  df.tallies <- df |> 
    select(Gene, x, afr_y, asj_y, eas_y, nfe_y, sas_y, mutability)
  
  # Calculate RVIS for each ancestry
  df <- df.tallies |> 
    mutate( 
      rvis_afr = CalcRVIS(df.tallies, "afr_y"),
      rvis_eas = CalcRVIS(df.tallies, "eas_y"),
      rvis_asj = CalcRVIS(df.tallies, "asj_y"),
      rvis_nfe = CalcRVIS(df.tallies, "nfe_y"),
      rvis_sas = CalcRVIS(df.tallies, "sas_y")
    )
  
  df_long <- df %>% 
    select(c(Gene, all_of(score_cols))) |> 
    pivot_longer(cols = score_cols, 
                 names_to = "score_type", 
                 values_to = "score_value")|>
    filter(!is.na(score_value))
  
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
