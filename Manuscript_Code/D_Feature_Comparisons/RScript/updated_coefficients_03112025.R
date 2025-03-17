library(ggplot2)
library(reshape2)
#install.packages("fillpattern")
library(fillpattern)
library(scales)
library(writexl)

#Top 20
non_zero_count_value = 20
output_path = '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/Top20/'

# #Top 1
# non_zero_count_value = 1
# output_path = '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/Top1/'

# #Top 5
# non_zero_count_value = 5
# output_path = '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/Top5/'

# #Top 10
# non_zero_count_value = 10
# output_path = '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/Top10/'

# #Top 50
# non_zero_count_value = 50
# output_path = '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/Top50/'

# #Top 3
# non_zero_count_value = 3
# output_path = '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/Top3/'

# Store datasets with fewer than 20 non-zero values dynamically
remove_names <- c()

# Modified process_files_by_pattern function
process_files_by_pattern <- function(pattern, directory = ".", feature_column = "Feature", importance_column = "Mean_Importance") {
  # Get a list of all files matching the pattern
  files <- list.files(path = directory, pattern = pattern, recursive = TRUE, full.names = TRUE)
  
  # Load all files into a list
  data_list <- lapply(files, read.csv)
  
  # Name the list elements with the corresponding filenames
  file_names <- basename(files)
  names(data_list) <- file_names
  
  # Process each dataframe
  processed_data_list <- lapply(names(data_list), function(file) {
    df <- data_list[[file]]
    if (all(c(feature_column, importance_column) %in% colnames(df))) {
      # Filter and sort the dataframe
      df <- df[order(-df[[importance_column]]), c(feature_column, importance_column), drop = FALSE]
      
      # Remove NA values from the importance column
      df <- df[!is.na(df[[importance_column]]), ]
      
      # Count non-zero values
      non_zero_count <- sum(df[[importance_column]] > 0, na.rm = TRUE)
      if (non_zero_count < non_zero_count_value) {
        datasets_with_few_nonzeros <<- c(datasets_with_few_nonzeros, basename(file))
      }
      
      return(df)
    } else {
      warning(paste("Required columns not found in:", file))
      return(NULL)
    }
  })
  
  # Assign file names to the list names
  processed_data_list <- setNames(processed_data_list, names(data_list))
  
  return(processed_data_list)
}

# Erawijantari
era_metab_results <- process_files_by_pattern("metab_importance\\.csv$", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Era_feature_extraction')
era_taxa_results <- process_files_by_pattern("taxa_importance\\.csv$", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Era_feature_extraction')
era_concat_results <- process_files_by_pattern("concat_importance\\.csv$", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Era_feature_extraction')
# Franzosa
franzosa_metab_results <- process_files_by_pattern("metab_importance\\.csv$", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Franzosa_feature_extraction')
franzosa_taxa_results <- process_files_by_pattern("taxa_importance.csv", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Franzosa_feature_extraction')
franzosa_concat_results <- process_files_by_pattern("concat_importance.csv", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Franzosa_feature_extraction')
#Yachida
yachida_metab_results <- process_files_by_pattern("metab_importance\\.csv$", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Yachida_feature_extraction')
yachida_taxa_results <- process_files_by_pattern("taxa_importance.csv", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Yachida_feature_extraction')
yachida_concat_results <- process_files_by_pattern("concat_importance.csv", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Yachida_feature_extraction')
#Wang
wang_metab_results <- process_files_by_pattern("metab_importance\\.csv$", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Wang_feature_extraction')
wang_taxa_results <- process_files_by_pattern("taxa_importance.csv", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Wang_feature_extraction')
wang_concat_results <- process_files_by_pattern("concat_importance.csv", '/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Coefficients/Wang_feature_extraction')



#split filtered and unfiltered dataframes
split_dataframes <- function(data_list, prefix_era, prefix_nf = "^nf") {
  # Ensure the list has names
  names_list <- names(data_list)
  # Split the data based on the prefixes
  filtered_df <- data_list[grepl(prefix_era, names_list)]
  unfiltered_df <- data_list[grepl(prefix_nf, names_list)]
  # Return a list with the two split datasets
  return(list(filtered = filtered_df, unfiltered = unfiltered_df))
}
#era
era_metab_split <- split_dataframes(era_metab_results, '^era', '^nf')
era_taxa_split <- split_dataframes(era_taxa_results, '^era', '^nf')
era_concat_split <- split_dataframes(era_concat_results, '^era', '^nf')
#franzosa
franz_metab_split <- split_dataframes(franzosa_metab_results, '^franzosa', '^un')
franz_taxa_split <- split_dataframes(franzosa_taxa_results, '^franzosa', '^un')
franz_concat_split <- split_dataframes(franzosa_concat_results, '^franzosa', '^un')
#yachida
yachida_metab_split <- split_dataframes(yachida_metab_results, '^yachida', '^nf')
yachida_taxa_split <- split_dataframes(yachida_taxa_results, '^yachida', '^nf')
yachida_concat_split <- split_dataframes(yachida_concat_results, '^yachida', '^nf')
#wang
wang_metab_split <- split_dataframes(wang_metab_results, '^wang', '^nf')
wang_taxa_split <- split_dataframes(wang_taxa_results, '^wang', '^nf')
wang_concat_split <- split_dataframes(wang_concat_results, '^wang', '^nf')

#metabolomics filtered and unfiltered
metab_df <- list(era_metab_split, franz_metab_split, yachida_metab_split, wang_metab_split)
metab_filtered_list <- do.call(c, lapply(metab_df, function(x) x$filtered))
metab_unfiltered_list <- do.call(c, lapply(metab_df, function(x) x$unfiltered))
#taxa unfiltered and filtered list
taxa_df <- list(era_taxa_split, franz_taxa_split, yachida_taxa_split, wang_taxa_split)
taxa_filtered_list <- do.call(c, lapply(taxa_df, function(x) x$filtered))
taxa_unfiltered_list <- do.call(c, lapply(taxa_df, function(x) x$unfiltered))
#concat unfiltered and filtered list
concat_df <- list(era_concat_split, franz_concat_split, yachida_concat_split, wang_concat_split)
concat_filtered_list <- do.call(c, lapply(concat_df, function(x) x$filtered))
concat_unfiltered_list <- do.call(c, lapply(concat_df, function(x) x$unfiltered))


# Print the names of datasets that will be removed
if (length(remove_names) > 0) {
  message("Datasets with fewer than Desired non-zero values in 'Mean_Importance' and will be removed:")
  print(remove_names)
}

# Filter the list to exclude the specified dataframes
metab_filtered_list1 <- metab_filtered_list[!(names(metab_filtered_list) %in% remove_names)]
metab_unfiltered_list1 <- metab_unfiltered_list[!(names(metab_unfiltered_list) %in% remove_names)]
taxa_filtered_list1 <- taxa_filtered_list[!(names(taxa_filtered_list) %in% remove_names)]
taxa_unfiltered_list1 <- taxa_unfiltered_list[!(names(taxa_unfiltered_list) %in% remove_names)]
concat_filtered_list1 <- concat_filtered_list[!(names(concat_filtered_list) %in% remove_names)]
concat_unfiltered_list1 <- concat_unfiltered_list[!(names(concat_unfiltered_list) %in% remove_names)]


# Function to extract top desired features from a list of dataframes for comparison
extract_top_features <- function(dataframes_list, feature_num) {
  top_features_list <- lapply(names(dataframes_list), function(name) {
    top_20 <- dataframes_list[[name]]$Feature[1:feature_num] # Extract top 20 features
  })
  # Name each list element with the original dataframe name
  names(top_features_list) <- names(dataframes_list)
  return(top_features_list)
}

metab_filtered_list2 <- extract_top_features(metab_filtered_list1, non_zero_count_value)
metab_unfiltered_list2 <- extract_top_features(metab_unfiltered_list1, non_zero_count_value)
taxa_filtered_list2 <- extract_top_features(taxa_filtered_list1, non_zero_count_value)
taxa_unfiltered_list2 <- extract_top_features(taxa_unfiltered_list1, non_zero_count_value)
concat_filtered_list2 <- extract_top_features(concat_filtered_list1, non_zero_count_value)
concat_unfiltered_list2 <- extract_top_features(concat_unfiltered_list1, non_zero_count_value)



# Convert list of lists to a dataframe with names as identifiers
list_to_dataframe <- function(results_list) {
  # Remove lists containing NA values
  valid_results <- results_list[!sapply(results_list, function(x) any(is.na(unlist(x))))]
  
  # Get names of valid lists
  valid_names <- names(valid_results)
  
  # Convert to dataframe
  df <- do.call(rbind, lapply(seq_along(valid_results), function(i) {
    cbind(Dataframe_Name = valid_names[i], as.data.frame(valid_results[[i]]))
  }))
  
  return(df)
}

separate_lists_by_pattern <- function(data_list, pattern) {
  # Filter the list based on the given pattern in the names
  matched_list <- data_list[grepl(pattern, names(data_list))]
  return(matched_list)
}
# Function to automate the process
automate_comparisons <- function(output_path, patterns, filtered_list, file_name) {
  
  # Ensure output_path exists
  #if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  # Separate lists based on patterns
  separated_lists <- lapply(patterns, function(p) separate_lists_by_pattern(filtered_list, paste0("^", p)))
  
  # Function to compute detailed comparisons
  compute_detailed_comparisons <- function(lists) {
    if (length(lists) != 3) return(NA)  # Skip if lists do not have exactly 3 elements
    
    enet <- lists[[1]]
    rf <- lists[[2]]
    xgb <- lists[[3]]
    
    format_count <- function(count) round(count / non_zero_count_value, 2)
    
    common_all <- intersect(intersect(enet, rf), xgb)
    enet_rf_overlap <- setdiff(intersect(enet, rf), xgb)
    enet_xgb_overlap <- setdiff(intersect(enet, xgb), rf)
    rf_xgb_overlap <- setdiff(intersect(rf, xgb), enet)
    
    enet_only <- setdiff(enet, union(rf, xgb))
    rf_only <- setdiff(rf, union(enet, xgb))
    xgb_only <- setdiff(xgb, union(enet, rf))
    
    list(
      all_overlap = format_count(length(common_all)),
      enet_rf_overlap = format_count(length(enet_rf_overlap)),
      enet_xgb_overlap = format_count(length(enet_xgb_overlap)),
      rf_xgb_overlap = format_count(length(rf_xgb_overlap)),
      enet_only = format_count(length(enet_only)),
      rf_only = format_count(length(rf_only)),
      xgb_only = format_count(length(xgb_only))
    )
  }
  # Apply the function to each separated list
  comparisons_results <- lapply(separated_lists, compute_detailed_comparisons)
  names(comparisons_results) <- patterns
  
  # Function to convert the result list to a dataframe
  convert_list_to_dataframe <- function(overlap_list) {
    data.frame(
      Method = c("all_overlap", "enet_rf", "enet_xgb", "rf_xgb", "enet", "rf", "xgb"),
      enet = c(overlap_list$all_overlap, overlap_list$enet_rf_overlap, overlap_list$enet_xgb_overlap, 0, overlap_list$enet_only, 0, 0),
      rf   = c(overlap_list$all_overlap, overlap_list$enet_rf_overlap, 0, overlap_list$rf_xgb_overlap, 0, overlap_list$rf_only, 0),
      xgb  = c(overlap_list$all_overlap, 0, overlap_list$enet_xgb_overlap, overlap_list$rf_xgb_overlap, 0, 0, overlap_list$xgb_only)
    )
  }
  for (name in names(comparisons_results)) {
    overlap_list <- comparisons_results[[name]]
    if (is.null(overlap_list) || all(is.na(overlap_list))) next  # Skip invalid results
    df <- convert_list_to_dataframe(overlap_list)
    #generate_and_save_plot(df, name)
  }
  #message("All plots have been successfully generated and saved in: ", output_path)
  comparisons_results <- list_to_dataframe(comparisons_results)
  writexl::write_xlsx(comparisons_results, path = paste0(output_path, non_zero_count_value, file_name, sep=''))
  return(comparisons_results)
}

extract_unique_patterns <- function(name_list) {
  # Define regex pattern to extract portion before numbers or model names
  pattern <- "^(.*?)(?:_\\d+|_enet|_xgb|_xgboost|_rf).*"
  # Extract the relevant part of the name
  extracted <- sub(pattern, "\\1", name_list)
  # Keep unique values
  unique_patterns <- unique(extracted)
  return(unique_patterns)
}

#pattern extraction for filtered metabolomics
fm_patterns <- extract_unique_patterns(names(metab_filtered_list2))
metab_filtered_results_list <- automate_comparisons(output_path,fm_patterns, metab_filtered_list2, "metab_filtered_03112025.xlsx")
#pattern extraction for unfiltered metabolomics
um_patterns <- extract_unique_patterns(names(metab_unfiltered_list2))
metab_unfiltered_results_list <- automate_comparisons(output_path, um_patterns, metab_unfiltered_list2, "metab_unfiltered_03112025.xlsx")
#pattern for filtered taxa
ft_patterns <- extract_unique_patterns(names(taxa_filtered_list2))
taxa_filtered_results_list <- automate_comparisons(output_path,ft_patterns, taxa_filtered_list2, "taxa_filtered_03112025.xlsx")
#pattern for unfiltered taxa
ut_patterns <- extract_unique_patterns(names(taxa_unfiltered_list2))
taxa_filtered_results_list <- automate_comparisons(output_path, ut_patterns, taxa_unfiltered_list2, "taxa_unfiltered_03112025.xlsx")
#pattern extraction for filtered concatenation
cf_patterns <- extract_unique_patterns(names(concat_filtered_list2))
concat_filtered_results_list <- automate_comparisons(output_path, cf_patterns, concat_filtered_list2, "concat_filtered_03112025.xlsx")
#pattern extraction for unfiltered concatenation
uc_patterns <- extract_unique_patterns(names(concat_unfiltered_list2))
concat_unfiltered_results_list <- automate_comparisons(output_path, uc_patterns, concat_unfiltered_list2, "concat_unfiltered_03112025.xlsx")

