library(dplyr)
library(tidyr)
# Function to read all CSV files from a directory and remove a specified column
read_all_csv <- function(directory_path, column_to_remove = NULL) {
  # List all CSV files in the specified directory
  files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Read all CSV files into a list of dataframes
  data_list <- lapply(files, function(file) {
    df <- read.csv(file)
    # Remove the specified column if it exists
    if (!is.null(column_to_remove) && column_to_remove %in% colnames(df)) {
      df <- df[, !colnames(df) %in% column_to_remove, drop = FALSE]
    }
    return(df)
  })
  names(data_list) <- basename(files)
  # Return the list of dataframes
  return(data_list)
}

# Load and clean data from various directories
franzosa_data <- read_all_csv(
  "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/franzosa_performance_metrics/results", 
  "Date")
wang_data <- read_all_csv(
  "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/wang_performance_csv/results", 
  "Date")
yachida_data <- read_all_csv(
  "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/yachida_performance_csv/results", 
  "Date")
era_data <- read_all_csv(
  "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Erawijantari_performance_csv/results", 
  "Date")


# Combine the lists into one
all_data <- c(franzosa_data, era_data, wang_data, yachida_data)

# Ensure the names of the list correspond to file names
names(all_data) <- basename(names(all_data))

# Filter dataframes whose names contain "filtered_testing_aucroc"

# Function to merge only common columns
merge_common_columns <- function(data_list) {
  # Identify common column names across all dataframes
  common_cols <- Reduce(intersect, lapply(data_list, colnames))
  
  # Select only the common columns from each dataframe
  filtered_data <- lapply(data_list, function(df) df[, common_cols, drop = FALSE])
  
  # Merge the filtered dataframes
  merged_df <- do.call(rbind, filtered_data)
  
  return(merged_df)
}

# Apply to filtered and unfiltered datasets
filtered_testing_aucroc <- merge_common_columns(all_data[grepl("_filtered_testing_aucroc", names(all_data))])
filtered_training_aucroc <- merge_common_columns(all_data[grepl("_filtered_training_aucroc", names(all_data))])
filtered_testing_rmse <- merge_common_columns(all_data[grepl("_filtered_testing_rmse", names(all_data))])
filtered_training_rmse <- merge_common_columns(all_data[grepl("_filtered_training_rmse", names(all_data))])
unfiltered_testing_aucroc <- merge_common_columns(all_data[grepl("unfiltered_testing_aucroc", names(all_data))])
unfiltered_training_aucroc <- merge_common_columns(all_data[grepl("unfiltered_training_aucroc", names(all_data))])
unfiltered_testing_rmse <- merge_common_columns(all_data[grepl("unfiltered_testing_rmse", names(all_data))])
unfiltered_training_rmse <- merge_common_columns(all_data[grepl("unfiltered_training_rmse", names(all_data))])

# Function to add a Dataset column based on specific keywords in File_Name
add_dataset_column <- function(df) {
  # Check for the presence of the File_Name column (case-insensitive)
  col_name <- names(df)[tolower(names(df)) == "file_name"]
  
  if (length(col_name) > 0 && nrow(df) > 0) {
    # Use the identified column name
    df$Dataset <- ifelse(grepl("franzosa", df[[col_name]], ignore.case = TRUE), "franzosa",
                         ifelse(grepl("yachida", df[[col_name]], ignore.case = TRUE), "yachida",
                                ifelse(grepl("wang", df[[col_name]], ignore.case = TRUE), "wang",
                                       ifelse(grepl("era", df[[col_name]], ignore.case = TRUE), "era", NA))))
  } else {
    warning("The 'File_Name' column is missing or the dataframe is empty.")
    df$Dataset <- NA  # Add a Dataset column with NA if File_Name is missing
  }
  return(df)
}

# Apply the function to your dataframes
filtered_testing_aucroc <- add_dataset_column(filtered_testing_aucroc)
filtered_training_aucroc <- add_dataset_column(filtered_training_aucroc)
filtered_testing_rmse <- add_dataset_column(filtered_testing_rmse)
filtered_training_rmse <- add_dataset_column(filtered_training_rmse)
unfiltered_testing_aucroc <- add_dataset_column(unfiltered_testing_aucroc)
unfiltered_training_aucroc <- add_dataset_column(unfiltered_training_aucroc)
unfiltered_testing_rmse <- add_dataset_column(unfiltered_testing_rmse)
unfiltered_training_rmse <- add_dataset_column(unfiltered_training_rmse)

convert_xgb_to_xgboost <- function(dataframe) {
  dataframe$Model <- gsub("^xgb$", "xgboost", dataframe$Model)
  return(dataframe)
}

# Applying the function to all specified dataframes
filtered_testing_aucroc <- convert_xgb_to_xgboost(filtered_testing_aucroc)
filtered_training_aucroc <- convert_xgb_to_xgboost(filtered_training_aucroc)
filtered_testing_rmse <- convert_xgb_to_xgboost(filtered_testing_rmse)
filtered_training_rmse <- convert_xgb_to_xgboost(filtered_training_rmse)
unfiltered_testing_aucroc <- convert_xgb_to_xgboost(unfiltered_testing_aucroc)
unfiltered_training_aucroc <- convert_xgb_to_xgboost(unfiltered_training_aucroc)
unfiltered_testing_rmse <- convert_xgb_to_xgboost(unfiltered_testing_rmse)
unfiltered_training_rmse <- convert_xgb_to_xgboost(unfiltered_training_rmse)


process_dataframe <- function(df) {
  # Retain the specified columns and combine Dataset, Response_Variable, and Model into one
  df <- df %>%
    select(Metabolomics, MSS, Concatenated, Averaged.Stacked, Weighted.NNLS, Lasso.Stacked, PLS, 
           Dataset, Response_Variable, Model) %>%
    mutate(Combined_Column = paste(Dataset, Response_Variable, Model, sep = "_")) %>%
    select(Metabolomics, MSS, Concatenated, Averaged.Stacked, Weighted.NNLS, Lasso.Stacked, PLS, Combined_Column)
  return(df)
}

# matrix created for matrix subtraction
filtered_testing_aucroc1 <- process_dataframe(filtered_testing_aucroc)
filtered_training_aucroc1 <- process_dataframe(filtered_training_aucroc)
filtered_testing_rmse1 <- process_dataframe(filtered_testing_rmse)
filtered_training_rmse1 <- process_dataframe(filtered_training_rmse)
unfiltered_testing_aucroc1 <- process_dataframe(unfiltered_testing_aucroc)
unfiltered_training_aucroc1 <- process_dataframe(unfiltered_training_aucroc)
unfiltered_testing_rmse1 <- process_dataframe(unfiltered_testing_rmse)
unfiltered_training_rmse1 <- process_dataframe(unfiltered_training_rmse)

# Define the function
update_combined_column <- function(data, column_name = "Combined_Column") {
  # Check if the specified column exists
  if (!column_name %in% colnames(data)) {
    stop(paste("Column", column_name, "does not exist in the dataframe."))
  }
  # Apply transformations
  data[[column_name]][1:6] <- paste0(data[[column_name]][1:6], "_ConvsCD")
  data[[column_name]][7:12] <- paste0(data[[column_name]][7:12], "_ConvsDisease")
  data[[column_name]][13:18] <- paste0(data[[column_name]][13:18], "_ConvsUC")
  # Return the modified dataframe
  return(data)
}

filtered_testing_rmse1 <- update_combined_column(filtered_testing_rmse1)
filtered_training_rmse1 <- update_combined_column(filtered_training_rmse1)
unfiltered_testing_rmse1 <- update_combined_column(unfiltered_testing_rmse1)
unfiltered_training_rmse1 <- update_combined_column(unfiltered_training_rmse1)

# Define the function
subtract_dataframes <- function(df_A, df_B) {
  # Ensure both dataframes have the same Combined_Column values
  all_combined_columns <- union(df_A$Combined_Column, df_B$Combined_Column)
  
  # Identify missing rows
  missing_in_A <- setdiff(all_combined_columns, df_A$Combined_Column)
  missing_in_B <- setdiff(all_combined_columns, df_B$Combined_Column)
  
  # Add missing rows with zeros to both dataframes
  df_A <- df_A %>%
    complete(Combined_Column = all_combined_columns, fill = list(Metabolomics = 0, MSS = 0, Concatenated = 0, 
                                                                 Averaged.Stacked = 0, Weighted.NNLS = 0, 
                                                                 Lasso.Stacked = 0, PLS = 0))
  
  df_B <- df_B %>%
    complete(Combined_Column = all_combined_columns, fill = list(Metabolomics = 0, MSS = 0, Concatenated = 0, 
                                                                 Averaged.Stacked = 0, Weighted.NNLS = 0, 
                                                                 Lasso.Stacked = 0, PLS = 0))
  
  # Sort both dataframes by Combined_Column to ensure alignment
  df_A <- df_A %>% arrange(Combined_Column)
  df_B <- df_B %>% arrange(Combined_Column)
  
  # Print rows where zeros were inserted
  if (length(missing_in_A) > 0) {
    cat("Rows added with zeros for df_A:\n")
    print(df_A %>% filter(Combined_Column %in% missing_in_A))
  }
  if (length(missing_in_B) > 0) {
    cat("Rows added with zeros for df_B:\n")
    print(df_B %>% filter(Combined_Column %in% missing_in_B))
  }
  
  # Perform element-wise subtraction on numeric columns only
  numeric_columns <- setdiff(names(df_A), "Combined_Column") # Exclude the "Combined_Column"
  result <- df_A
  result[, numeric_columns] <- df_A[, numeric_columns] - df_B[, numeric_columns]
  
  # Return the resulting dataframe
  return(result)
}

testing_aucroc_difference <- subtract_dataframes(filtered_testing_aucroc1, unfiltered_testing_aucroc1)
testing_rmse_difference <- subtract_dataframes(unfiltered_testing_rmse1, filtered_testing_rmse1)
training_aucroc_difference <- subtract_dataframes(filtered_training_aucroc1, unfiltered_training_aucroc1)
training_rmse_difference <- subtract_dataframes(unfiltered_testing_rmse1, filtered_testing_rmse1)

# Load required library
library(writexl)

# Example usage
write_xlsx(list(testing_aucroc_difference = testing_aucroc_difference, 
                testing_rmse_difference = testing_rmse_difference,
                training_aucroc_difference = training_aucroc_difference,
                training_rmse_difference = training_rmse_difference), "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/performance_difference_12242024.xlsx")


########Heatmap

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)

# Function to create a heatmap
create_heatmap <- function(data, heatmap_name, metric = "AUCROC") {
  # Define numeric columns to include in the heatmap
  numeric_columns <- c("Metabolomics", "MSS", "Concatenated", 
                       "Averaged.Stacked", "Weighted.NNLS", 
                       "Lasso.Stacked", "PLS")
  
  # Subset data to include only numeric columns
  heatmap_data <- data[, numeric_columns]
  
  # Prepare annotations
  row_annotation <- data.frame(
    Response_Variable = data$Response_Variable,
    Model = data$Model,
    Dataset = data$Dataset
  )
  
  # Ensure unique row names for ComplexHeatmap
  row.names(heatmap_data) <- paste(data$Response_Variable, seq_len(nrow(data)), sep = "_")
  
  # Annotation for colors
  model_colors <- c("enet" = "white", "rf" = "gray", "xgboost" = "black")
  dataset_colors <- c("franzosa" = "yellow", "wang" = "blue", "yachida" = "purple", "era" = "green")
  
  # Set color scale based on metric
  if (metric == "AUCROC") {
    # AUCROC: White (0) to Red (1)
    custom_colors <- colorRamp2(c(0, 1), c("white", "red"))
    heatmap_data_scaled <- heatmap_data  # No scaling needed for AUCROC
  } else if (metric == "RMSE") {
    # RMSE: Log10 scale for color mapping
    max_value <- max(heatmap_data, na.rm = TRUE)  # Calculate max RMSE value
    min_positive <- min(heatmap_data[heatmap_data > 0], na.rm = TRUE)  # Minimum positive value
    custom_colors <- colorRamp2(
      c(log10(min_positive), log10(max_value)),  # Include all positive values in log10 range
      c("white", "red")
    )
    # Transform data for color mapping (but keep original values for display)
    heatmap_data_scaled <- log10(pmax(heatmap_data, min_positive))  # Avoid log10(0)
  } else {
    stop("Unsupported metric. Use 'AUCROC' or 'RMSE'.")
  }
  
  # Create row annotation
  row_ha <- rowAnnotation(
    Model = row_annotation$Model,
    Dataset = row_annotation$Dataset,
    col = list(
      Model = model_colors,
      Dataset = dataset_colors
    )
  )
  
  # Create the heatmap
  Heatmap(
    as.matrix(if (metric == "RMSE") heatmap_data_scaled else heatmap_data),
    name = heatmap_name,  # Use provided heatmap name
    col = custom_colors, 
    row_split = row_annotation$Dataset,  # Split rows by Dataset
    row_title = NULL,  # Suppress row titles
    cluster_rows = FALSE,  # Do not cluster rows
    cluster_columns = FALSE,  # Do not cluster columns
    show_row_names = TRUE,  # Display row names
    row_names_gp = gpar(fontsize = 10),  # Adjust row name font size
    row_labels = row_annotation$Response_Variable,  # Use Response_Variable for y-axis labels
    left_annotation = row_ha,  # Add the annotations
    border = TRUE,   # Add white borders between cells
    rect_gp = gpar(col = "white", lwd = 2),
    cell_fun = function(j, i, x, y, width, height, fill) {
      # Display original values in cells
      grid.text(sprintf("%.2f", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
    }
  )
}

# Define the default path for saving PDFs
output_path <- "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/heatmap_figure2"

# Save high-quality heatmaps as PDFs
pdf(file.path(output_path, "Filtered_Testing_AUCROC.pdf"), width = 8, height = 8)
create_heatmap(filtered_testing_aucroc, "Filtered Testing AUCROC", metric = "AUCROC")
dev.off()

pdf(file.path(output_path, "Filtered_Training_AUCROC.pdf"), width = 8, height = 8)
create_heatmap(filtered_training_aucroc, "Filtered Training AUCROC", metric = "AUCROC")
dev.off()

pdf(file.path(output_path, "Unfiltered_Testing_AUCROC.pdf"), width = 8, height = 8)
create_heatmap(unfiltered_testing_aucroc, "Unfiltered Testing AUCROC", metric = "AUCROC")
dev.off()

pdf(file.path(output_path, "Unfiltered_Training_AUCROC.pdf"), width = 8, height = 8)
create_heatmap(unfiltered_training_aucroc, "Unfiltered Training AUCROC", metric = "AUCROC")
dev.off()

pdf(file.path(output_path, "Filtered_Testing_RMSE.pdf"), width = 8, height = 8)
create_heatmap(filtered_testing_rmse, "log10(Filtered Testing RMSE)", metric = "RMSE")
dev.off()

pdf(file.path(output_path, "Filtered_Training_RMSE.pdf"), width = 8, height = 8)
create_heatmap(filtered_training_rmse, "log10(Filtered Training RMSE)", metric = "RMSE")
dev.off()

pdf(file.path(output_path, "Unfiltered_Testing_RMSE.pdf"), width = 8, height = 6)
create_heatmap(unfiltered_testing_rmse, "log10(Unfiltered Testing RMSE)", metric = "RMSE")
dev.off()

pdf(file.path(output_path, "Unfiltered_Training_RMSE.pdf"), width = 8, height = 6)
create_heatmap(unfiltered_training_rmse, "log10(Unfiltered Training RMSE)", metric = "RMSE")
dev.off()


