# Function to normalize importance by absolute values
normalize_absolute <- function(df) {
  df_normalized <- as.data.frame(lapply(df, function(x) {
    abs(x) / max(abs(x), na.rm = TRUE)
  }))
  return(df_normalized)
}

# ---------------------------------------------------------------
# Choose the correct feature-importance routine at run-time
# ---------------------------------------------------------------
if (tolower(model_to_run) == "elastic net") {
  
  feature_extraction <- function(model_list) {
    importance_list <- list()
    
    for (i in seq_along(model_list)) {
      current_model   <- model_list[[i]]
      final_glmnet    <- current_model$finalModel
      optimal_lambda  <- final_glmnet$lambdaOpt
      coefficients    <- coef(final_glmnet, s = optimal_lambda)
      
      feature_df <- as.data.frame(as.matrix(coefficients))
      feature_df <- data.frame(
        Feature   = rownames(feature_df),
        Importance = feature_df[, 1]
      )
      rownames(feature_df) <- NULL
      feature_df <- feature_df[feature_df$Feature != "(Intercept)", ]
      colnames(feature_df)[2] <- paste0("Importance_Model_", i)
      
      importance_list[[paste0("Model_", i)]] <- feature_df
    }
    
    combined <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE),
                       importance_list)
    normalized <- normalize_absolute(combined[, -1])
    out <- cbind(Feature = combined$Feature, normalized)
    out$Mean_Importance <- rowMeans(out[, -1], na.rm = TRUE)
    out
  }
  
} else if (tolower(model_to_run) == "random forest") {
  
  feature_extraction <- function(model_list) {
    importance_list <- list()
    
    for (i in seq_along(model_list)) {
      current_model <- model_list[[i]]
      feature_df    <- as.data.frame(current_model$finalModel$importance)
      feature_df$Feature <- rownames(feature_df)
      rownames(feature_df) <- NULL
      colnames(feature_df)[1] <- paste0("Importance_Model_", i)
      
      importance_list[[paste0("Model_", i)]] <- feature_df
    }
    
    combined   <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE),
                         importance_list)
    normalized <- normalize_absolute(combined[, -1])
    out <- cbind(Feature = combined$Feature, normalized)
    out$Mean_Importance <- rowMeans(out[, -1], na.rm = TRUE)
    out
  }
  
} else if (tolower(model_to_run) == "xgboost") {
  
  feature_extraction <- function(model_list) {
    importance_list <- list()
    
    for (i in seq_along(model_list)) {
      current_model <- model_list[[i]]
      final_model   <- current_model$finalModel
      
      feature_df <- xgb.importance(
        model = final_model,
        feature_names = colnames(current_model$trainingData)[
          -ncol(current_model$trainingData)
        ]
      )
      feature_df <- feature_df[, c("Feature", "Gain")]
      colnames(feature_df)[2] <- paste0("Importance_Model_", i)
      
      importance_list[[paste0("Model_", i)]] <- feature_df
    }
    
    combined            <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE),
                                  importance_list)
    combined[is.na(combined)] <- 0
    normalized          <- normalize_absolute(combined[, -1])
    out <- cbind(Feature = combined$Feature, normalized)
    out$Mean_Importance <- rowMeans(out[, -1], na.rm = TRUE)
    out
  }
  
} else {
  stop("Unknown model type. Valid options are 'elastic net', 'random forest', or 'xgboost'.")
}