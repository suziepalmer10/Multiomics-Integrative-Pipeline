# Load required packages
library(tools)  # For file path manipulation
library(caret)
library(xgboost)

# Function to normalize importance by absolute values
normalize_absolute <- function(df) {
  df_normalized <- as.data.frame(lapply(df, function(x) {
    abs(x) / max(abs(x), na.rm = TRUE)
  }))
  return(df_normalized)
}

calculate_normalized_importance_enet <- function(model_list) {
  importance_list <- list()
  
  for (i in seq_along(model_list)) {
    current_model <- model_list[[i]]
    final_glmnet <- current_model$finalModel
    optimal_lambda <- final_glmnet$lambdaOpt
    coefficients <- coef(final_glmnet, s = optimal_lambda)
    feature_importance <- as.data.frame(as.matrix(coefficients))
    feature_importance <- data.frame(Feature = rownames(feature_importance), Importance = feature_importance[, 1])
    rownames(feature_importance) <- NULL
    feature_importance <- feature_importance[feature_importance$Feature != "(Intercept)", ]
    colnames(feature_importance)[2] <- paste0("Importance_Model_", i)
    importance_list[[paste0("Model_", i)]] <- feature_importance
  }
  
  combined_importance <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE), importance_list)
  importance_only <- combined_importance[, -1]
  normalized_importance <- normalize_absolute(importance_only)
  combined_importance_abs_normalized <- cbind(Feature = combined_importance$Feature, normalized_importance)
  combined_importance_abs_normalized$Mean_Importance <- rowMeans(combined_importance_abs_normalized[, -1], na.rm = TRUE)
  return(combined_importance_abs_normalized)
}

calculate_normalized_importance_rf <- function(model_list) {
  # Initialize an empty list to store the feature importance for each model
  importance_list <- list()
  # Loop through each model in the provided model list
  for (i in seq_along(model_list)) {
    # Extract the current model
    current_model <- model_list[[i]]
    # Use varImp to get the feature importance from the random forest model
    feature_importance <- current_model$finalModel$importance
    # Convert to a data frame for easier handling
    feature_importance <- as.data.frame(feature_importance)
    feature_importance$Feature <- rownames(feature_importance)
    rownames(feature_importance) <- NULL
    # Rename the importance column to represent the current model
    colnames(feature_importance)[1] <- paste0("Importance_Model_", i)
    # Store in the list with model name as identifier
    importance_list[[paste0("Model_", i)]] <- feature_importance
  }
  # Combine all feature importance data frames into one using 'Feature' column
  combined_importance <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE), importance_list)
  # Extract importance values and normalize
  importance_only <- combined_importance[, -1]
  normalized_importance <- normalize_absolute(importance_only)
  # Add the 'Feature' column back to the normalized dataframe
  combined_importance_abs_normalized <- cbind(Feature = combined_importance$Feature, normalized_importance)
  # Calculate the mean importance across all models
  combined_importance_abs_normalized$Mean_Importance <- rowMeans(combined_importance_abs_normalized[, -1], na.rm = TRUE)
  # Return the normalized and combined feature importance dataframe
  return(combined_importance_abs_normalized)
}

calculate_normalized_importance_xgb <- function(model_list) {
  # Initialize an empty list to store the feature importance for each model
  importance_list <- list()
  # Loop through each model in the provided model list
  for (i in seq_along(model_list)) {
    #print(i)
    # Extract the current model
    current_model <- model_list[[i]]
    final_model <- current_model$finalModel
    # Extract feature importance using xgb.importance for xgbTree
    feature_importance <- xgb.importance(model = final_model, feature_names = colnames(current_model$trainingData)[-ncol(current_model$trainingData)])
    feature_importance <- as.data.frame(feature_importance[, c("Feature", "Gain")])
    #print(nrow(feature_importance))
    colnames(feature_importance)[2] <- "Importance"
    # Rename the importance column to represent the current model
    colnames(feature_importance)[which(names(feature_importance) != "Feature")] <- paste0("Importance_Model_", i)
    # Store in the list with model name as identifier
    importance_list[[paste0("Model_", i)]] <- feature_importance
  }
  # Combine all feature importance data frames into one using 'Feature' column
  combined_importance <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE), importance_list)
  combined_importance[is.na(combined_importance)] <- 0
  # Extract importance values and normalize
  importance_only <- combined_importance[, -1]
  normalized_importance <- normalize_absolute(importance_only)
  # Add the 'Feature' column back to the normalized dataframe
  combined_importance_abs_normalized <- cbind(Feature = combined_importance$Feature, normalized_importance)
  # Calculate the mean importance across all models
  combined_importance_abs_normalized$Mean_Importance <- rowMeans(combined_importance_abs_normalized[, -1], na.rm = TRUE)
  # Return the normalized and combined feature importance dataframe
  return(combined_importance_abs_normalized)
}

create_statistics_df <- function(df, column_name) {
  breaks <- seq(0, 1, by = 0.1)
  bins <- cut(df[[column_name]], breaks, include.lowest = TRUE)
  statistics <- as.data.frame(table(bins))
  count_zero <- sum(df[[column_name]] == 0)
  new_row <- data.frame(bins = "0", Freq = count_zero)
  statistics <- rbind(new_row, statistics)
  return(statistics)
}

process_rdata_file <- function(file_path, model_kind) {
  # Construct the full file path if not absolute
  full_path <- normalizePath(file_path, mustWork = TRUE)
  # Extract base name without extension
  base_name <- gsub("\\.RData$", "", basename(full_path))
  print(base_name)
  # Load the RData file
  load(full_path)
  # Process models (assumes metab_model, concat_model, and mss_model exist in RData file)
  if (model_kind == 'enet') {
    metab_importance <- calculate_normalized_importance_enet(metab_model)
    concat_importance <- calculate_normalized_importance_enet(concat_model)
    taxa_importance <- calculate_normalized_importance_enet(mss_model)
  } else if (model_kind == 'rf') {
    metab_importance <- calculate_normalized_importance_rf(metab_model)
    concat_importance <- calculate_normalized_importance_rf(concat_model)
    taxa_importance <- calculate_normalized_importance_rf(mss_model)
  } else if (model_kind == 'xgb') {
    metab_importance <- calculate_normalized_importance_xgb(metab_model)
    concat_importance <- calculate_normalized_importance_xgb(concat_model)
    taxa_importance <- calculate_normalized_importance_xgb(mss_model)
  }
  
  metab_statistics <- create_statistics_df(metab_importance, "Mean_Importance")
  concat_statistics <- create_statistics_df(concat_importance, "Mean_Importance")
  taxa_statistics <- create_statistics_df(taxa_importance, "Mean_Importance")
  
  # Write output files
  write.csv(metab_importance, file = paste0(base_name, "_metab_importance.csv"), row.names = FALSE)
  write.csv(concat_importance, file = paste0(base_name, "_concat_importance.csv"), row.names = FALSE)
  write.csv(taxa_importance, file = paste0(base_name, "_taxa_importance.csv"), row.names = FALSE)
  write.csv(metab_statistics, file = paste0(base_name, "_metab_statistics.csv"), row.names = FALSE)
  write.csv(concat_statistics, file = paste0(base_name, "_concat_statistics.csv"), row.names = FALSE)
  write.csv(taxa_statistics, file = paste0(base_name, "_taxa_statistics.csv"), row.names = FALSE)
}

extract_importance <- function(input_dir, pattern_, model_kind_) {
  # Get all files ending with ".RData"
  rdata_files <- list.files(path = input_dir, pattern = pattern_, full.names = TRUE)
  print(length(rdata_files))
  # Process each file
  for (file_path in rdata_files) {
    process_rdata_file(file_path, model_kind_)
  }
}
#franzosa_filtered
setwd('/work/PCDC/s180020/IntegratedLearner_Dec2024/Franzosa_integrated/IntegratedLearner/results/RDataFiles/Franzosa_filtered_RData')
input_dir <- "/work/PCDC/s180020/IntegratedLearner_Dec2024/Franzosa_integrated/IntegratedLearner/results/RDataFiles/Franzosa_filtered_RData"
extract_importance(input_dir, "enet.RData", 'enet')
extract_importance(input_dir, "xgboost.RData", 'xgb')
extract_importance(input_dir, "rf.RData", 'rf')
#franzosa_unfiltered
setwd('/work/PCDC/s180020/IntegratedLearner_Dec2024/Franzosa_integrated/IntegratedLearner/results/RDataFiles/franzosa_original')
input_dir <- "/work/PCDC/s180020/IntegratedLearner_Dec2024/Franzosa_integrated/IntegratedLearner/results/RDataFiles/franzosa_original"
extract_importance(input_dir, "enet.RData", 'enet')
extract_importance(input_dir, "xgboost.RData", 'xgb')
extract_importance(input_dir, "rf.RData", 'rf')
#wang_filtered_and_unfiltered
setwd('/work/PCDC/s180020/IntegratedLearner_Dec2024/Wang_integrated/IntegratedLearner/results/RDataFiles')
input_dir <- "/work/PCDC/s180020/IntegratedLearner_Dec2024/Wang_integrated/IntegratedLearner/results/RDataFiles"
extract_importance(input_dir, "enet.RData", 'enet')
extract_importance(input_dir, "xgb.RData", 'xgb')
extract_importance(input_dir, "rf.RData", 'rf')
#era filtered and unfiltered
setwd('/work/PCDC/s180020/IntegratedLearner_Dec2024/Erawijantari_integrated/IntegratedLearner/results/RDataFiles')
input_dir <- "/work/PCDC/s180020/IntegratedLearner_Dec2024/Erawijantari_integrated/IntegratedLearner/results/RDataFiles"
extract_importance(input_dir, "enet.RData", 'enet')
extract_importance(input_dir, "xgb.RData", 'xgb')
extract_importance(input_dir, "rf.RData", 'rf')
#yachida filtered and unfiltered
setwd('/work/PCDC/s180020/IntegratedLearner_Dec2024/Yachida_integrated/IntegratedLearner/results/RDataFiles')
input_dir <- "/work/PCDC/s180020/IntegratedLearner_Dec2024/Yachida_integrated/IntegratedLearner/results/RDataFiles"
extract_importance(input_dir, "enet.RData", 'enet')
extract_importance(input_dir, "xgb.RData", 'xgb')
extract_importance(input_dir, "rf.RData", 'rf')