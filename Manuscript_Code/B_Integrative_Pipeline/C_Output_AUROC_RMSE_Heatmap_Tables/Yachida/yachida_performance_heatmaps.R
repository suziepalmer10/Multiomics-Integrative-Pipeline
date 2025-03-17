#This code is to extract performance metrics to be used for heatmap analysis
#Use for unfiltered and filtered franzosa data.
# Load necessary libraries
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(grid)

# Set the directory containing your CSV files

csv_directory <- "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/yachida_performance_csv"
# Load in CSV files
csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)

# Create a dataframe to store metadata extracted from filenames
metadata <- data.frame(File_Name = basename(csv_files), stringsAsFactors = FALSE)

# Process files starting with "wang"
filtered_df <- metadata %>%
  filter(str_starts(File_Name, "yachida")) %>%
  mutate(
    Study_Name = str_split_fixed(File_Name, "_", 4)[, 1],
    Response_Variable = str_split_fixed(File_Name, "_", 4)[, 2],
    Date = str_split_fixed(File_Name, "_", 4)[, 3],
    Model = str_replace(str_split_fixed(File_Name, "_", 4)[, 4], ".csv", "")
  )

# Process files starting with "nf"
unfiltered_df <- metadata %>%
  filter(str_starts(File_Name, "nf")) %>%
  mutate(
    Study_Name = str_split_fixed(File_Name, "_", 5)[, 2],
    Response_Variable = str_split_fixed(File_Name, "_", 5)[, 3],
    Date = str_split_fixed(File_Name, "_", 5)[, 4],
    Model = str_replace(str_split_fixed(File_Name, "_", 5)[, 5], ".csv", "")
  )

# Function to read CSV files based on filtered dataframe
read_csv_files <- function(metadata_df, csv_directory) {
  # Get the full file paths from the metadata
  file_paths <- file.path(csv_directory, metadata_df$File_Name)
  
  # Read each file into a list of dataframes
  file_contents <- lapply(file_paths, read.csv, stringsAsFactors = FALSE)
  
  # Add metadata to each dataframe as an attribute
  names(file_contents) <- metadata_df$File_Name
  
  return(file_contents)
}

filtered_files <- read_csv_files(filtered_df, csv_directory)
unfiltered_files <- read_csv_files(unfiltered_df, csv_directory)


extract_training_metrics <- function(csv_data_list, metric_columns, result_columns) {
  # Initialize an empty dataframe to store extracted values
  extracted_values <- data.frame(matrix(ncol = length(result_columns), nrow = length(csv_data_list)))
  colnames(extracted_values) <- result_columns
  
  # Loop through each dataframe in the list
  for (i in seq_along(csv_data_list)) {
    data <- csv_data_list[[i]]
    
    # Check if the primary metric column exists
    if (metric_columns[1] %in% colnames(data)) {
      # Extract the metric column and transpose
      extracted_values[i, ] <- t(data[[metric_columns[1]]][1:length(result_columns)])
    } 
    # Check for the alternative metric column
    else if (metric_columns[2] %in% colnames(data)) {
      extracted_values[i, ] <- t(data[[metric_columns[2]]][1:length(result_columns)])
    } else {
      # Fill with NA if neither column exists
      extracted_values[i, ] <- NA
    }
  }
  
  return(extracted_values)
}

# Define the RMSE and AUC-ROC column names to search for
training_metric_columns <- c("RMSE.Train._Mean", "Aucroc.Train._Mean")
testing_metric_columns <- c("RMSE.Test._Mean", "Aucroc.Test._Mean")
# Define the result dataframe's column names
rmse_columns <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked", 
                  "Weighted NNLS", "Lasso Stacked", "PLS")

unique(filtered_df$Response_Variable)
continuous_var <- c("age", "alcohol", "BMI", "BrinkmanIndex")
binary_var <- c("gender", "healthyvscancer", "healthyvsearly", "healthyvsstageI.II", "healthyvsstageIII.IV")

#Filtered Data
filtered_training <- extract_training_metrics(filtered_files, training_metric_columns, rmse_columns)
filtered_training <- cbind(filtered_df, filtered_training)
filtered_testing <- extract_training_metrics(filtered_files, testing_metric_columns, rmse_columns)
filtered_testing <- cbind(filtered_df, filtered_testing)
unique(filtered_training$Response_Variable)
filtered_training <- filtered_training %>% arrange(Study_Name, Response_Variable, Model)
filtered_testing <- filtered_testing %>% arrange(Study_Name, Response_Variable, Model)
#subset by continuous or binary
filtered_training_aucroc <- subset(filtered_training, Response_Variable %in% binary_var)
filtered_training_rmse  <- subset(filtered_training, Response_Variable %in% continuous_var)
filtered_testing_aucroc <- subset(filtered_testing, Response_Variable %in% binary_var)
filtered_testing_rmse  <- subset(filtered_testing, Response_Variable %in% continuous_var)
#Unfiltered Data
unfiltered_training <- extract_training_metrics(unfiltered_files, training_metric_columns, rmse_columns)
unfiltered_training <- cbind(unfiltered_df, unfiltered_training)
unfiltered_testing <- extract_training_metrics(unfiltered_files, testing_metric_columns, rmse_columns)
unfiltered_testing <- cbind(unfiltered_df, unfiltered_testing)
unique(unfiltered_training$Response_Variable)
unfiltered_training <- unfiltered_training %>% arrange(Study_Name, Response_Variable, Model)
unfiltered_testing <- unfiltered_testing %>% arrange(Study_Name, Response_Variable, Model)
#subset by continuous or binary
unfiltered_training_aucroc <- subset(unfiltered_training, Response_Variable %in% binary_var)
unfiltered_training_rmse  <- subset(unfiltered_training, Response_Variable %in% continuous_var)
unfiltered_testing_aucroc <- subset(unfiltered_testing, Response_Variable %in% binary_var)
unfiltered_testing_rmse  <- subset(unfiltered_testing, Response_Variable %in% continuous_var)

PATH <- "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/yachida_performance_csv/results"

write.csv(unfiltered_training_aucroc, file.path(PATH, "yachida_unfiltered_training_aucroc_12152024.csv"), row.names = FALSE)
write.csv(filtered_training_aucroc, file.path(PATH,"yachida_filtered_training_aucroc_12152024.csv"), row.names = FALSE)
write.csv(unfiltered_testing_aucroc, file.path(PATH,"yachida_unfiltered_testing_aucroc_12152024.csv"), row.names = FALSE)
write.csv(filtered_testing_aucroc, file.path(PATH,"yachida_filtered_testing_aucroc_12152024.csv"), row.names = FALSE)

write.csv(unfiltered_training_rmse, file.path(PATH,"yachida_unfiltered_training_rmse_12152024.csv"), row.names = FALSE)
write.csv(filtered_training_rmse, file.path(PATH,"yachida_filtered_training_rmse_12152024.csv"), row.names = FALSE)
write.csv(unfiltered_testing_rmse, file.path(PATH,"yachida_unfiltered_testing_rmse_12152024.csv"), row.names = FALSE)
write.csv(filtered_testing_rmse, file.path(PATH,"yachida_filtered_testing_rmse_12152024.csv"), row.names = FALSE)




# Function to create the heatmap
create_heatmap <- function(df, heatmap_name, title) {
  # Combine the columns into a new y-axis label
  df <- df %>%
    mutate(Y_Axis_Labels = paste(Study_Name, Response_Variable, Model, sep = "_"))

  
  # Prepare the data for the heatmap
  heatmap_data <- df[, rmse_columns]
  heatmap_data <- as.matrix(heatmap_data)  # Ensure it's a matrix
  rownames(heatmap_data) <- df$Y_Axis_Labels
  
  # Define the colors for the "Model" column
  model_colors <- c("enet" = "white", "rf" = "gray", "xgb" = "black")
  
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

# Function to save each heatmap as an individual PDF
save_heatmap_as_pdf <- function(data, metric, title) {
  # Create a valid filename from the title
  filename <- file.path(output_directory, paste0(gsub(" ", "_", title), ".pdf")) 
  #filename <- paste0(gsub(" ", "_", title), ".pdf")  # Replace spaces with underscores for the filename
  # Open a PDF device with the generated filename
  pdf(filename, width = 8, height = 6)  # Adjust dimensions as needed
  # Create the heatmap
  create_heatmap(data, metric, title)
  # Close the PDF device
  dev.off()
}

output_directory = "/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/yachida_performance_csv"
# Call the function for each heatmap
save_heatmap_as_pdf(filtered_training_aucroc, 'AUC-ROC Mean', "Filtered AUC-ROC Training")
save_heatmap_as_pdf(unfiltered_training_aucroc, 'AUC-ROC Mean', "Unfiltered AUC-ROC Training")
save_heatmap_as_pdf(filtered_testing_aucroc, 'AUC-ROC Mean', "Filtered AUC-ROC Testing")
save_heatmap_as_pdf(unfiltered_testing_aucroc, 'AUC-ROC Mean', "Unfiltered AUC-ROC Testing")
save_heatmap_as_pdf(filtered_training_rmse, 'RMSE Mean', "Filtered RMSE Training")
save_heatmap_as_pdf(unfiltered_training_rmse, 'RMSE Mean', "Unfiltered RMSE Training")
save_heatmap_as_pdf(filtered_testing_rmse, 'RMSE Mean', "Filtered RMSE Testing")
save_heatmap_as_pdf(unfiltered_testing_rmse, 'RMSE Mean', "Unfiltered RMSE Testing")





