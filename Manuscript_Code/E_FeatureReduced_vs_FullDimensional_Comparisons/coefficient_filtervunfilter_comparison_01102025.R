library(writexl)

process_csv_files <- function(path, pattern, comparison_grep, top_features) {
  # Get a list of all CSV files matching the pattern
  csv_files <- list.files(path, pattern = pattern, full.names = TRUE)
  
  # Load all CSV files into a named list
  data_list <- setNames(lapply(csv_files, function(file) {
    data <- read.csv(file)
    return(data)
  }), tools::file_path_sans_ext(basename(csv_files)))
  
  # Keep only "Feature" and "Mean_Importance" columns
  data_list <- lapply(data_list, function(df) {
    df <- df[, c("Feature", "Mean_Importance"), drop = FALSE]
    return(df)
  })
  
  # Keep only the top N elements based on Mean_Importance
  data_list <- lapply(data_list, function(df) {
    df <- df[order(-df$Mean_Importance), ]  # Order by descending Mean_Importance
    df <- head(df, top_features)            # Keep top N rows
    return(df)
  })
  
  # Check and remove dataframes with NULL values or zero Mean_Importance
  data_list1 <- lapply(names(data_list), function(name) {
    df <- data_list[[name]]
    if (any(is.na(df$Mean_Importance))) {
      cat(name, "has been removed because it contains NA values in Mean_Importance.\n")
      return(NULL)
    }
    if (any(df$Mean_Importance == 0)) {
      cat(name, "has been removed because it contains zero Mean_Importance values.\n")
      return(NULL)
    }
    return(df)
  })
  names(data_list1) <- names(data_list)
  data_list1 <- Filter(Negate(is.null), data_list1)
  
  matches <- list()
  unmatched <- c()
  
  for (name in names(data_list1)) {
    clean_name <- sub("\\..*$", "", name)
    
    if (startsWith(clean_name, "nf_")) {
      match_name <- sub("^nf_", "", clean_name)
    } else if (startsWith(clean_name, "un_")) {
      match_name <- sub("^un_", "", clean_name)
    } else {
      match_name_nf <- paste0("nf_", clean_name)
      match_name_un <- paste0("un_", clean_name)
    }
    
    possible_matches <- sub("\\..*$", "", names(data_list1))
    if (exists("match_name") && match_name %in% possible_matches) {
      matches[[clean_name]] <- match_name
    } else if (exists("match_name_nf") && match_name_nf %in% possible_matches) {
      matches[[clean_name]] <- match_name_nf
    } else if (exists("match_name_un") && match_name_un %in% possible_matches) {
      matches[[clean_name]] <- match_name_un
    } else {
      unmatched <- c(unmatched, clean_name)
    }
  }
  
  # Compare Features for matching pairs and calculate the ratio
  comparison_results <- sapply(names(matches), function(clean_name) {
    df1 <- data_list1[[clean_name]]
    df2 <- data_list1[[matches[[clean_name]]]]
    
    common_features <- intersect(df1$Feature, df2$Feature)
    
    comparison_name <- paste0(clean_name, "_vs_nf")
    
    ratio <- length(common_features) / top_features
    
    return(setNames(ratio, comparison_name))
  })
  
  comparison_results <- unlist(comparison_results)
  
  if (length(unmatched) > 0) {
    cat("The following dataframes do not have a match:\n")
    print(unmatched)
  } else {
    cat("All dataframes have matches.\n")
  }
  
  comparison_df <- data.frame(
    Dataset = names(comparison_results),
    Proportion = as.numeric(comparison_results),
    stringsAsFactors = FALSE
  )
  
  comparison_df$Dataset <- sub("^.*\\.", "", comparison_df$Dataset)
  comparison_df <- comparison_df[!grepl(comparison_grep, comparison_df$Dataset), ]
  
  enet_df <- comparison_df[grep("enet", comparison_df$Dataset), ]
  rf_df <- comparison_df[grep("rf", comparison_df$Dataset), ]
  xgb_df <- comparison_df[grep("xgb", comparison_df$Dataset), ]
  
  return(list(enet_df = enet_df, rf_df = rf_df, xgb_df = xgb_df))
}

setwd('/Users/suzettepalmer/Desktop/Integrated_Pipeline_March2025/Code_and_Analyses/Manuscript_Code/C_Feature_Importance_Values/Output_FeatureExtraction/')

feature_thresholds <- c(1, 5, 10, 20, 50)

base_dir <- "/Users/suzettepalmer/Desktop/Integrated_Pipeline_March2025/Code_and_Analyses/Manuscript_Code/E_FeatureReduced_vs_FullDimensional_Comparisons/Comparisons/"

for (top_n in feature_thresholds) {
  #dir_name <- paste0("Top", top_n)
  dir_name <- file.path(base_dir, paste0("Top", top_n))
  dir.create(dir_name, showWarnings = FALSE)
  
  # Process metabolite, taxa, and concatenated results
  era_metabolite_result <- process_csv_files("Era_feature_extraction/", "metab_importance.csv$", "^nf", top_n)
  era_taxa_result <- process_csv_files("Era_feature_extraction", "taxa_importance.csv$", "^nf", top_n)
  era_concat_result <- process_csv_files("Era_feature_extraction", "concat_importance.csv$", "^nf", top_n)
  
  franz_metabolite_result <- process_csv_files("Franzosa_feature_extraction/", "metab_importance.csv$", "^un", top_n)
  franz_taxa_result <- process_csv_files("Franzosa_feature_extraction/", "taxa_importance.csv$", "^un", top_n)
  franz_concat_result <- process_csv_files("Franzosa_feature_extraction/", "concat_importance.csv$", "^un", top_n)
  
  wang_metabolite_result <- process_csv_files("Wang_feature_extraction/", "metab_importance.csv$", "^nf", top_n)
  wang_taxa_result <- process_csv_files("Wang_feature_extraction/", "taxa_importance.csv$", "^nf", top_n)
  wang_concat_result <- process_csv_files("Wang_feature_extraction/", "concat_importance.csv$", "^nf", top_n)
  
  yachida_metabolite_result <- process_csv_files("Yachida_feature_extraction/", "metab_importance.csv$", "^nf", top_n)
  yachida_taxa_result <- process_csv_files("Yachida_feature_extraction/", "taxa_importance.csv$", "^nf", top_n)
  yachida_concat_result <- process_csv_files("Yachida_feature_extraction/", "concat_importance.csv$", "^nf", top_n)
  
  # Create result lists for metabolites, taxa, and concatenated datasets
  metabolite_results <- list(
    "Elastic Net" = do.call(rbind, lapply(list(
      era_metabolite_result$enet_df, franz_metabolite_result$enet_df,
      wang_metabolite_result$enet_df, yachida_metabolite_result$enet_df
    ), as.data.frame)),
    
    "Random Forest" = do.call(rbind, lapply(list(
      era_metabolite_result$rf_df, franz_metabolite_result$rf_df,
      wang_metabolite_result$rf_df, yachida_metabolite_result$rf_df
    ), as.data.frame)),
    
    "XGBoost" = do.call(rbind, lapply(list(
      era_metabolite_result$xgb_df, franz_metabolite_result$xgb_df,
      wang_metabolite_result$xgb_df, yachida_metabolite_result$xgb_df
    ), as.data.frame))
  )
  
  taxa_results <- list(
    "Elastic Net" = do.call(rbind, lapply(list(
      era_taxa_result$enet_df, franz_taxa_result$enet_df,
      wang_taxa_result$enet_df, yachida_taxa_result$enet_df
    ), as.data.frame)),
    
    "Random Forest" = do.call(rbind, lapply(list(
      era_taxa_result$rf_df, franz_taxa_result$rf_df,
      wang_taxa_result$rf_df, yachida_taxa_result$rf_df
    ), as.data.frame)),
    
    "XGBoost" = do.call(rbind, lapply(list(
      era_taxa_result$xgb_df, franz_taxa_result$xgb_df,
      wang_taxa_result$xgb_df, yachida_taxa_result$xgb_df
    ), as.data.frame))
  )
  
  concat_results <- list(
    "Elastic Net" = do.call(rbind, lapply(list(
      era_concat_result$enet_df, franz_concat_result$enet_df,
      wang_concat_result$enet_df, yachida_concat_result$enet_df
    ), as.data.frame)),
    
    "Random Forest" = do.call(rbind, lapply(list(
      era_concat_result$rf_df, franz_concat_result$rf_df,
      wang_concat_result$rf_df, yachida_concat_result$rf_df
    ), as.data.frame)),
    
    "XGBoost" = do.call(rbind, lapply(list(
      era_concat_result$xgb_df, franz_concat_result$xgb_df,
      wang_concat_result$xgb_df, yachida_concat_result$xgb_df
    ), as.data.frame))
  )
  
  write_xlsx(metabolite_results, file.path(dir_name, paste0("metabolite_comp_results_Top", top_n, ".xlsx")))
  write_xlsx(taxa_results, file.path(dir_name, paste0("taxa_comp_results_Top", top_n, ".xlsx")))
  write_xlsx(concat_results, file.path(dir_name, paste0("concat_comp_results_Top", top_n, ".xlsx")))

}
