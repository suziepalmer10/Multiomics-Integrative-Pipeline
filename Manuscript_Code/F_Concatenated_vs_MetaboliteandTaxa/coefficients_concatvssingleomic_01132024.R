library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(writexl)


#top_n = 20

setwd('/Users/suzettepalmer/Desktop/Integrated_Pipeline_March2025/Code_and_Analyses/Manuscript_Code/C_Feature_Importance_Values/Output_FeatureExtraction')

feature_thresholds <- c(1, 5, 10, 20, 50)

base_dir <- "/Users/suzettepalmer/Desktop/Integrated_Pipeline_March2025/Code_and_Analyses/Manuscript_Code/F_Concatenated_vs_MetaboliteandTaxa"

for (top_n in feature_thresholds) {
  dir_name <- file.path(base_dir, paste0("Top", top_n))
  dir.create(dir_name, showWarnings = FALSE)
# Global list to store datasets with fewer than 20 features
removed_datasets <- list()

# Function to process CSV files
process_csv_files <- function(path, pattern) {
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
  
  # Keep only the top n elements based on Mean_Importance
  data_list <- lapply(data_list, function(df) {
    df <- df[order(-df$Mean_Importance), ]  # Order by descending Mean_Importance
    df <- head(df, top_n)                      # Keep top n rows
    return(df)
  })
  
  # Check and remove dataframes with NULL values, zero Mean_Importance, or fewer than 20 features
  data_list1 <- lapply(names(data_list), function(name) {
    df <- data_list[[name]]
    
    # Remove datasets with NA values in Mean_Importance
    if (any(is.na(df$Mean_Importance))) {
      cat(name, "has been removed because it contains NA values in Mean_Importance.\n")
      removed_datasets[[name]] <<- df  # Store in removed datasets list
      return(NULL)
    }
    
    # Remove datasets with Mean_Importance values of zero
    if (any(df$Mean_Importance == 0)) {
      cat(name, "has been removed because it contains zero Mean_Importance values.\n")
      removed_datasets[[name]] <<- df  # Store in removed datasets list
      return(NULL)
    }
    
    # Remove datasets with fewer than n features
    if (nrow(df) < top_n) {
      cat(name, "has been removed because it has fewer than n features.\n")
      removed_datasets[[name]] <<- df  # Store in removed datasets list
      return(NULL)
    }
    
    return(df)
  })
  
  names(data_list1) <- names(data_list)
  
  # Remove NULL dataframes from the list
  data_list1 <- Filter(Negate(is.null), data_list1)
  
  return(data_list1)
}

# Process each dataset
# Erawijantari
era_metabolite_result <- process_csv_files("Era_feature_extraction/", "metab_importance.csv$")
era_taxa_result <- process_csv_files("Era_feature_extraction/", "taxa_importance.csv$")
era_concat_result <- process_csv_files("Era_feature_extraction/", "concat_importance.csv$")

# Franzosa
franz_metabolite_result <- process_csv_files("Franzosa_feature_extraction/", "metab_importance.csv$")
franz_taxa_result <- process_csv_files("Franzosa_feature_extraction/", "taxa_importance.csv$")
franz_concat_result <- process_csv_files("Franzosa_feature_extraction/", "concat_importance.csv$")

# Wang
wang_metabolite_result <- process_csv_files("Wang_feature_extraction/", "metab_importance.csv$")
wang_taxa_result <- process_csv_files("Wang_feature_extraction/", "taxa_importance.csv$")
wang_concat_result <- process_csv_files("Wang_feature_extraction/", "concat_importance.csv$")

# Yachida
yachida_metabolite_result <- process_csv_files("Yachida_feature_extraction/", "metab_importance.csv$")
yachida_taxa_result <- process_csv_files("Yachida_feature_extraction/", "taxa_importance.csv$")
yachida_concat_result <- process_csv_files("Yachida_feature_extraction/", "concat_importance.csv$")

# Print datasets removed due to < n features
if (length(removed_datasets) > 0) {
  message("\nDatasets removed due to having fewer than n features:")
  print(names(removed_datasets))
}

#This counts the number of metabolites and taxa in the concatenated files. 
# Function to count features starting with "t__" and "m__"
count_features <- function(df) {
  taxa_count <- sum(grepl("^t__", df$Feature))
  metabolites_count <- sum(grepl("^m__", df$Feature))
  return(data.frame(Taxa = taxa_count, Metabolites = metabolites_count))
}

# Apply the function to each dataframe in the list and combine the results into one dataframe
concat_era_count <- do.call(rbind, lapply(era_concat_result, count_features))
concat_franzosa_count <- do.call(rbind, lapply(franz_concat_result, count_features))
concat_wang_count <- do.call(rbind, lapply(wang_concat_result, count_features))
concat_yachida_count <- do.call(rbind, lapply(yachida_concat_result, count_features))

preprocess_for_bar_plot <- function(df, taxa_col = "Taxa", metabolite_col = "Metabolites") {
  # Add row_id as a new column
  df$row_id <- rownames(df)
  
  # Create the 'Model' column based on patterns in 'row_id'
  df$Model <- ifelse(grepl("enet", df$row_id), "enet",
                     ifelse(grepl("rf", df$row_id), "rf",
                            ifelse(grepl("xgb", df$row_id), "xgb", "other")))
  return(df)
  
}

preprocess_era_count <- preprocess_for_bar_plot(concat_era_count)
preprocess_franz_count <- preprocess_for_bar_plot(concat_franzosa_count)
preprocess_wang_count <- preprocess_for_bar_plot(concat_wang_count)
preprocess_yachida_count <- preprocess_for_bar_plot(concat_yachida_count)

#rbind the datasets
preprocess_datasets <- rbind(preprocess_era_count, preprocess_franz_count, 
                             preprocess_wang_count, preprocess_yachida_count)

#subset filtered and unfiltered datasets.
subset_unfiltered <- preprocess_datasets[grepl("^(nf|un)", preprocess_datasets$row_id), ]
subset_filtered <- preprocess_datasets[!grepl("^(nf|un)", preprocess_datasets$row_id), ]

# Function to export datasets into an Excel file with separate sheets for each model
export_to_excel <- function(data, filename) {
  # Split data into three lists based on the 'Model' column
  model_splits <- split(data, data$Model)
  # Keep only relevant models
  model_splits <- model_splits[c("enet", "rf", "xgb")]
  # Write the list to an Excel file
  #write_xlsx(model_splits, path = filename)
  write_xlsx(model_splits, file.path(dir_name, paste0(top_n, filename, ".xlsx")))
  cat("Exported:", filename, "\n")
}

# Export subset_unfiltered to an Excel file
export_to_excel(subset_unfiltered, "fulldimensional_concat_proportions")

# Export subset_filtered to an Excel file
export_to_excel(subset_filtered, "featurereduced_concat_proportions")

#Now for overlapping feature analysis
find_overlapping_features <- function(metab_list, concat_list, suffix) {
  metab_prefixes <- sub(suffix, "", names(metab_list))
  concat_prefixes <- sub("_concat_importance", "", names(concat_list))
  # Identify matching prefixes
  common_prefixes <- intersect(metab_prefixes, concat_prefixes)
  # Initialize a dataframe to store overlapping counts
  overlapping_features_df <- data.frame(Prefix = character(), Overlap_Count = numeric())
  # Loop through common prefixes to find overlaps
  for (prefix in common_prefixes) {
    metab_df <- metab_list[[paste0(prefix, suffix)]]
    concat_df <- concat_list[[paste0(prefix, "_concat_importance")]]
    
    # Find overlapping features in the "Feature" column
    overlap <- intersect(metab_df$Feature, concat_df$Feature)
    
    # Calculate count divided by n
    overlap_count <- length(overlap) / top_n
    
    # Append result to the dataframe
    overlapping_features_df <- rbind(overlapping_features_df, data.frame(Prefix = prefix, Overlap_Count = overlap_count))
  }
  return(overlapping_features_df)
}

#Era
era_overlaps_metab <- find_overlapping_features(era_metabolite_result, era_concat_result, "_metab_importance")
era_overlaps_taxa<- find_overlapping_features(era_taxa_result, era_concat_result, "_taxa_importance")
merged_era_overlaps <- merge(era_overlaps_metab, era_overlaps_taxa, by = "Prefix", all = TRUE)
merged_era_overlaps$Model <- gsub(".*_(.*)$", "\\1", merged_era_overlaps$Prefix)
merged_era_overlaps$Dataset <- gsub("_(?!.*_).*$", "", merged_era_overlaps$Prefix, perl = TRUE)
colnames(merged_era_overlaps) <- c("Prefix", "Metab_Concat Overlap", "Taxa_Concat Overlap", "Model", "Dataset")
nf_merged_era_overlaps <- merged_era_overlaps[grepl("^nf_", merged_era_overlaps$Prefix), ]
f_merged_era_overlaps <- merged_era_overlaps[!grepl("^nf_", merged_era_overlaps$Prefix), ]

#franzosa
franz_overlaps_metab <- find_overlapping_features(franz_metabolite_result, franz_concat_result, "_metab_importance")
franz_overlaps_taxa<- find_overlapping_features(franz_taxa_result, franz_concat_result, "_taxa_importance")
merged_franz_overlaps <- merge(franz_overlaps_metab, franz_overlaps_taxa, by = "Prefix", all = TRUE)
merged_franz_overlaps$Model <- gsub(".*_(.*)$", "\\1", merged_franz_overlaps$Prefix)
merged_franz_overlaps$Dataset <- gsub("_(?!.*_).*$", "", merged_franz_overlaps$Prefix, perl = TRUE)
colnames(merged_franz_overlaps) <- c("Prefix", "Metab_Concat Overlap", "Taxa_Concat Overlap", "Model", "Dataset")
merged_franz_overlaps$Dataset <- gsub("Convs", "Conv", merged_franz_overlaps$Dataset)
nf_merged_franz_overlaps <- merged_franz_overlaps[grepl("^un_", merged_franz_overlaps$Prefix), ]
f_merged_franz_overlaps <- merged_franz_overlaps[!grepl("^un_", merged_franz_overlaps$Prefix), ]

#wang
wang_overlaps_metab <- find_overlapping_features(wang_metabolite_result, wang_concat_result, "_metab_importance")
wang_overlaps_taxa<- find_overlapping_features(wang_taxa_result, wang_concat_result, "_taxa_importance")
merged_wang_overlaps <- merge(wang_overlaps_metab, wang_overlaps_taxa, by = "Prefix", all = TRUE)
merged_wang_overlaps$Model <- gsub(".*_(.*)$", "\\1", merged_wang_overlaps$Prefix)
merged_wang_overlaps$Dataset <- gsub("_(?!.*_).*$", "", merged_wang_overlaps$Prefix, perl = TRUE)
colnames(merged_wang_overlaps) <- c("Prefix", "Metab_Concat Overlap", "Taxa_Concat Overlap", "Model", "Dataset")
nf_merged_wang_overlaps <- merged_wang_overlaps[grepl("^nf_", merged_wang_overlaps$Prefix), ]
f_merged_wang_overlaps <- merged_wang_overlaps[!grepl("^nf_", merged_wang_overlaps$Prefix), ]

#yachida
yachida_overlaps_metab <- find_overlapping_features(yachida_metabolite_result, yachida_concat_result, "_metab_importance")
yachida_overlaps_taxa<- find_overlapping_features(yachida_taxa_result, yachida_concat_result, "_taxa_importance")
merged_yachida_overlaps <- merge(yachida_overlaps_metab, yachida_overlaps_taxa, by = "Prefix", all = TRUE)
merged_yachida_overlaps$Model <- gsub(".*_(.*)$", "\\1", merged_yachida_overlaps$Prefix)
merged_yachida_overlaps$Dataset <- gsub("_(?!.*_).*$", "", merged_yachida_overlaps$Prefix, perl = TRUE)
colnames(merged_yachida_overlaps) <- c("Prefix", "Metab_Concat Overlap", "Taxa_Concat Overlap", "Model", "Dataset")
nf_merged_yachida_overlaps <- merged_yachida_overlaps[grepl("^nf_", merged_yachida_overlaps$Prefix), ]
f_merged_yachida_overlaps <- merged_yachida_overlaps[!grepl("^nf_", merged_yachida_overlaps$Prefix), ]

#merge all datasets together
f_merged <- rbind(f_merged_era_overlaps, f_merged_franz_overlaps, f_merged_wang_overlaps, f_merged_yachida_overlaps)
nf_merged <- rbind(nf_merged_era_overlaps, nf_merged_franz_overlaps, nf_merged_wang_overlaps, nf_merged_yachida_overlaps)
#replace xgboost with xgb
f_merged <- f_merged %>% mutate(Model = ifelse(Model == "xgboost", "xgb", Model))
nf_merged <- nf_merged %>% mutate(Model = ifelse(Model == "xgboost", "xgb", Model))

# Export subset_unfiltered to an Excel file
export_to_excel(f_merged, "feature_reduced_comparisons")

# Export subset_filtered to an Excel file
export_to_excel(nf_merged, "full_dimensional_comparisons")
}

