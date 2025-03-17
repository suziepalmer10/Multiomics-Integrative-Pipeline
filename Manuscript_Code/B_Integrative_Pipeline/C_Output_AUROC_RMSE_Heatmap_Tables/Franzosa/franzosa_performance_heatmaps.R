#This code is to extract performance metrics to be used for heatmap analysis
#Use for unfiltered and filtered franzosa data.
# Load necessary libraries
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(grid)

# Set the directory containing your CSV files
#csv_directory <- "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/franzosa_performance_metrics/franzosa_original"
csv_directory <- "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/franzosa_performance_metrics/franzosa_filtered"

# List all CSV files in the directory
csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)

# Create a dataframe to store metadata extracted from filenames
metadata <- data.frame(File_Name = basename(csv_files))

# Separate the filename into different components
metadata <- metadata %>%
  mutate(
    Study_name = str_split_fixed(File_Name, "_", 5)[, 1],
    Subset = str_split_fixed(File_Name, "_", 5)[, 2],
    Response_Variable = str_split_fixed(File_Name, "_", 5)[, 3],
    Date = str_split_fixed(File_Name, "_", 5)[, 4],
    Model = str_replace(str_split_fixed(File_Name, "_", 5)[, 5], ".csv", "")
  )

# Function to read CSV files
read_csv_files <- function(files) {
  lapply(files, read.csv, stringsAsFactors = FALSE)
}


# Read in all CSV files into a list
csv_data_list <- read_csv_files(csv_files)

# Column names for RMSE Mean values
rmse_columns <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked", 
                  "Weighted NNLS", "Lasso Stacked", "PLS")

#Extract Training Values for RMSE and AUCROC
# Initialize an empty dataframe to store RMSE Mean values
rmse_values <- data.frame(matrix(ncol = length(rmse_columns), nrow = length(csv_data_list)))
colnames(rmse_values) <- rmse_columns
# Extract RMSE Mean or AUC-ROC from each CSV and add it to the df. 
for (i in seq_along(csv_data_list)) {
  # Assuming csv_data_list[[i]] is a dataframe and RMSE Mean is a column
  data <- csv_data_list[[i]]
  
  # Check if "RMSE Mean" column exists
  if ("RMSE.Train._Mean" %in% colnames(data)) {
    # Extract RMSE Mean column and transpose
    rmse_values[i, ] <- t(data[["RMSE.Train._Mean"]][1:length(rmse_columns)])
  } 
  else if ("Aucroc.Train._Mean"%in% colnames(data)) {
    # Fill with NA if the column doesn't exist
    rmse_values[i, ] <- t(data[["Aucroc.Train._Mean"]][1:length(rmse_columns)])
  }
}
final_metadata <- cbind(metadata, rmse_values)
unique(final_metadata$Response_Variable)
final_metadata$Response_Variable <- gsub("ConvCD", "ConvsCD", final_metadata$Response_Variable)
final_metadata$Response_Variable <- gsub("ConvUC", "ConvsUC", final_metadata$Response_Variable)
final_metadata$Response_Variable <- gsub("ConvDisease", "ConvsDisease", final_metadata$Response_Variable)
final_metadata <- final_metadata %>%
  arrange(Subset, Response_Variable, Model)

unfiltered_binary_training_aucroc <- subset(final_metadata, Response_Variable %in% c("ConvsCD", "ConvsUC", "ConvsDisease"))
unfiltered_continuous_training_rmse  <- subset(final_metadata, Response_Variable %in% c("Fp", "Age"))

#Extract TestingValues for RMSE and AUCROC
# Initialize an empty dataframe to store RMSE Mean values
rmse_values_test <- data.frame(matrix(ncol = length(rmse_columns), nrow = length(csv_data_list)))
colnames(rmse_values_test) <- rmse_columns
# Extract RMSE Mean or AUC-ROC from each CSV and add it to the df. 
for (i in seq_along(csv_data_list)) {
  # Assuming csv_data_list[[i]] is a dataframe and RMSE Mean is a column
  data <- csv_data_list[[i]]
  
  # Check if "RMSE Mean" column exists
  if ("RMSE.Test._Mean" %in% colnames(data)) {
    # Extract RMSE Mean column and transpose
    rmse_values_test[i, ] <- t(data[["RMSE.Test._Mean"]][1:length(rmse_columns)])
  } 
  else if ("Aucroc.Test._Mean"%in% colnames(data)) {
    # Fill with NA if the column doesn't exist
    rmse_values_test[i, ] <- t(data[["Aucroc.Test._Mean"]][1:length(rmse_columns)])
  }
}
final_metadata_test <- cbind(metadata, rmse_values_test)
unique(final_metadata_test$Response_Variable)
final_metadata_test$Response_Variable <- gsub("ConvCD", "ConvsCD", final_metadata_test$Response_Variable)
final_metadata_test$Response_Variable <- gsub("ConvUC", "ConvsUC", final_metadata_test$Response_Variable)
final_metadata_test$Response_Variable <- gsub("ConvDisease", "ConvsDisease", final_metadata_test$Response_Variable)
final_metadata_test <- final_metadata_test %>%
  arrange(Subset, Response_Variable, Model)

unfiltered_binary_test_aucroc <- subset(final_metadata_test, Response_Variable %in% c("ConvsCD", "ConvsUC", "ConvsDisease"))
unfiltered_continuous_test_rmse  <- subset(final_metadata_test, Response_Variable %in% c("Fp", "Age"))

# Define the base path
PATH <- "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/franzosa_performance_metrics/results"

# # Write the dataframes to the specified path
# write.csv(unfiltered_binary_test_aucroc, file.path(PATH, "franzosa_unfiltered_testing_aucroc_12152024.csv"), row.names = FALSE)
# write.csv(unfiltered_binary_training_aucroc, file.path(PATH, "franzosa_unfiltered_training_aucroc_12152024.csv"), row.names = FALSE)
# write.csv(unfiltered_continuous_test_rmse, file.path(PATH, "franzosa_unfiltered_testing_rmse_12152024.csv"), row.names = FALSE)
# write.csv(unfiltered_continuous_training_rmse, file.path(PATH, "franzosa_unfiltered_training_rmse_12152024.csv"), row.names = FALSE)
# Write the dataframes to the specified path
write.csv(unfiltered_binary_test_aucroc, file.path(PATH, "franzosa_filtered_testing_aucroc_12152024.csv"), row.names = FALSE)
write.csv(unfiltered_binary_training_aucroc, file.path(PATH, "franzosa_filtered_training_aucroc_12152024.csv"), row.names = FALSE)
write.csv(unfiltered_continuous_test_rmse, file.path(PATH, "franzosa_filtered_testing_rmse_12152024.csv"), row.names = FALSE)
write.csv(unfiltered_continuous_training_rmse, file.path(PATH, "franzosa_filtered_training_rmse_12152024.csv"), row.names = FALSE)




# Function to create the heatmap
create_heatmap <- function(df, heatmap_name, title) {
  # Combine the columns into a new y-axis label
  df <- df %>%
    mutate(Y_Axis_Labels = paste(Study_name, Subset, Response_Variable, Model, sep = "_"))

  
  # Prepare the data for the heatmap
  heatmap_data <- df[, rmse_columns]
  heatmap_data <- as.matrix(heatmap_data)  # Ensure it's a matrix
  rownames(heatmap_data) <- df$Y_Axis_Labels
  
  # Define the colors for the "Model" column
  model_colors <- c("enet" = "white", "rf" = "gray", "xgboost" = "black")
  
  # Create a row annotation for the "Model" column
  row_annotation <- rowAnnotation(
    Model = df$Model,  # Use the Model column as data
    annotation_legend_param = list(title = "Model"),  # Add a legend for clarity
    col = list(Model = model_colors)  # Map the colors to the Model values
  )
  
  # Create and return the heatmap
  ht <- Heatmap(
    heatmap_data,
    name = heatmap_name,
    row_labels = rownames(heatmap_data), # Use the combined labels as row labels
    show_row_names = TRUE,
    show_column_names = TRUE,
    #column_title = "Models",
    cluster_rows = FALSE,  # Optional: Disable clustering of rows
    cluster_columns = TRUE,  # Optional: Cluster columns
    left_annotation = row_annotation,  # Add the row annotation
    cell_fun = function(j, i, x, y, width, height, fill) {
      # Draw the filled rectangle with white borders
      grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white", lwd = 1))
      # Add text to the center of each cell with two decimal places
      grid.text(sprintf("%.2f", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
    }
  )
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = title)
}




#unfiltered
# create_heatmap(unfiltered_binary_test_aucroc, 'AUC-ROC Mean', "Unfiltered AUC-ROC Testing")
# create_heatmap(unfiltered_binary_training_aucroc, 'AUC-ROC Mean', "Unfiltered AUC-ROC Training")
# create_heatmap(unfiltered_continuous_test_rmse, 'RMSE Mean', "Unfiltered RMSE Testing")
# create_heatmap(unfiltered_continuous_training_rmse, 'RMSE Mean', "Unfiltered RMSE Training")
# create_heatmap(unfiltered_binary_test_aucroc, 'AUC-ROC Mean', "Filtered AUC-ROC Testing")
# create_heatmap(unfiltered_binary_training_aucroc, 'AUC-ROC Mean', "Filtered AUC-ROC Training")
# create_heatmap(unfiltered_continuous_test_rmse, 'RMSE Mean', "Filtered RMSE Testing")
# create_heatmap(unfiltered_continuous_training_rmse, 'RMSE Mean', "Filtered RMSE Training")



# Optionally, save the final dataframe as a CSV
write.csv(final_metadata, file = "final_metadata_with_rmse.csv", row.names = FALSE)

