library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

library(stringr)

# Define the base directory containing the "TopX" folders
base_dir <- "/Users/suzettepalmer/Desktop/Integrated_Pipeline_March2025/Code_and_Analyses/Manuscript_Code/E_FeatureReduced_vs_FullDimensional_Comparisons/Comparisons/"

# Feature thresholds (directories to iterate through)
feature_thresholds <- c(1, 5, 10, 20, 50)

# Define the lists for continuous and binary variables
continuous <- list(
  "era_age", "era_cholesterol", "era_glucose", 
  "franzosa_ControlvsCD_Age",  "franzosa_ControlvsCD_Fp",
  "franzosa_ControlvsDisease_Age", "franzosa_ControlvsDisease_Fp",
  "franzosa_ControlvsUC_Age", "franzosa_ControlvsUC_Fp",
  "yachida_age", "yachida_alcohol", "yachida_BrinkmanIndex", 
  "wang_age", "wang_creatinine", "wang_egfr", "wang_urea",
  "wang_bmi", "yachida_BMI"
)

binary <- list(
  "era_alcohol", "era_sg", "era_gender",
  "franzosa_ControlvsCD_ConvCD", 
  "franzosa_ControlvsDisease_ConvDisease",
  "franzosa_ControlvsUC_ConvUC",
  "yachida_gender",
  "yachida_healthyvscancer", "yachida_healthyvsearly", "yachida_healthyvsstageI_II", "yachida_healthyvsstageIII_IV",
  "wang_studygroup", "wang_gender"
)

# Initialize empty dataframes for each dataset type
metabolomics_data <- data.frame()
taxa_data <- data.frame()
concatenated_data <- data.frame()

# Loop through each TopX directory
for (top_n in feature_thresholds) {
  dir_name <- file.path(base_dir, paste0("Top", top_n))
  
  # Define paths to Excel files
  metabolite_file <- file.path(dir_name, paste0("metabolite_comp_results_Top", top_n, ".xlsx"))
  taxa_file <- file.path(dir_name, paste0("taxa_comp_results_Top", top_n, ".xlsx"))
  concat_file <- file.path(dir_name, paste0("concat_comp_results_Top", top_n, ".xlsx"))
  
  # List of files to process
  file_list <- list(
    Metabolomics = metabolite_file,
    Taxa = taxa_file,
    Concatenated = concat_file
  )
  
  # Iterate through each file type
  for (data_type in names(file_list)) {
    file_path <- file_list[[data_type]]
    
    if (file.exists(file_path)) {
      # Read each sheet separately
      elastic_net <- read_xlsx(file_path, sheet = "Elastic Net") %>%
        mutate(Model = "Elastic Net")
      rf <- read_xlsx(file_path, sheet = "Random Forest") %>%
        mutate(Model = "Random Forest")
      xgb <- read_xlsx(file_path, sheet = "XGBoost") %>%
        mutate(Model = "XGBoost")
      
      # Combine all sheets into one dataframe
      file_data <- bind_rows(elastic_net, rf, xgb) %>%
        mutate(Feature_Threshold = paste0("Top", top_n),
               Data_Type = data_type)
      
      # Store in the appropriate dataframe
      if (data_type == "Metabolomics") {
        metabolomics_data <- bind_rows(metabolomics_data, file_data)
      } else if (data_type == "Taxa") {
        taxa_data <- bind_rows(taxa_data, file_data)
      } else if (data_type == "Concatenated") {
        concatenated_data <- bind_rows(concatenated_data, file_data)
      }
    }
  }
}

# Ensure correct column names
colnames(metabolomics_data) <- c("Dataset", "Proportion", "Model", "Feature_Threshold", "Data_Type")
colnames(taxa_data) <- c("Dataset", "Proportion", "Model", "Feature_Threshold", "Data_Type")
colnames(concatenated_data) <- c("Dataset", "Proportion", "Model", "Feature_Threshold", "Data_Type")

# Convert columns to factors for correct ordering
for (df in list(metabolomics_data, taxa_data, concatenated_data)) {
  df$Feature_Threshold <- factor(df$Feature_Threshold, levels = paste0("Top", feature_thresholds))
  df$Model <- factor(df$Model, levels = c("Elastic Net", "Random Forest", "XGBoost"))
}


# Function to check if a dataset matches any value in a list
matches_list <- function(dataset, list_values) {
  sapply(dataset, function(x) any(str_detect(x, str_c(list_values, collapse = "|"))))
}

# Subset data using partial string matching
metabolomics_binary <- metabolomics_data[matches_list(metabolomics_data$Dataset, binary), ]
metabolomics_continuous <- metabolomics_data[matches_list(metabolomics_data$Dataset, continuous), ]

taxa_binary <- taxa_data[matches_list(taxa_data$Dataset, binary), ]
taxa_continuous <- taxa_data[matches_list(taxa_data$Dataset, continuous), ]

concatenated_binary <- concatenated_data[matches_list(concatenated_data$Dataset, binary), ]
concatenated_continuous <- concatenated_data[matches_list(concatenated_data$Dataset, continuous), ]

# Identify missing variables that were not categorized
all_datasets <- unique(c(metabolomics_data$Dataset, taxa_data$Dataset, concatenated_data$Dataset))
classified_datasets <- unique(c(metabolomics_binary$Dataset, metabolomics_continuous$Dataset,
                                taxa_binary$Dataset, taxa_continuous$Dataset,
                                concatenated_binary$Dataset, concatenated_continuous$Dataset))
unclassified_datasets <- setdiff(all_datasets, classified_datasets)




custom_colors <- c("Elastic Net" = "white",  # Blue
                   "Random Forest" = "#939598",  # Orange
                   "XGBoost" = "#B31E2D")  # Green

# Function to create and save violin plots
save_violin_plot <- function(data, filename, title) {
  data$Model <- factor(data$Model, levels = c("Elastic Net", "Random Forest", "XGBoost"))
  data$Feature_Threshold <- factor(data$Feature_Threshold, levels = c("Top1", "Top5", "Top10", "Top20", "Top50"))
  if (nrow(data) > 0) {  # Only create a plot if there is data
    p <- ggplot(data, aes(x = Model, y = Proportion, fill = Model)) +
      
      #geom_violin(alpha = 0.6) +
      #geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +  # Add points for individual values
      geom_violin(alpha = 0.7, trim = TRUE, width = 0.6) +
      geom_jitter(aes(color = Metric), position = position_jitter(width = 0.2, height = 0), 
                  alpha = 0.6, size = 1.5, shape = 1, color = 'black') + 
      facet_wrap(~Feature_Threshold, ncol = 1, scales = "fixed") +  # Facet into rows by TopX category
      #labs(title = title, x = "Model", y = "Proportion of Overlapping Features") +
      labs(title = title, 
           x = "Metric", y = "Proportion of Features") +
      
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.05, "lines"),
      ) + 
      #theme(legend.position = "none",
      #      strip.text = element_text(size = 12, face = "bold"))
      scale_fill_manual(values = custom_colors) +  # Apply custom colors for boxplot
      scale_y_continuous(limits = c(0, 1.0), expand = expansion(mult = c(0.075, 0.075))) 
    ggsave(filename, plot = p, width = 4, height = 6, dpi = 300)
  }
}


# Define output directory for plots
plot_dir <- file.path(base_dir, "Violin_Plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save the six violin plots as PDFs
save_violin_plot(metabolomics_binary, file.path(plot_dir, "Metabolomics_Binary_ViolinPlot.pdf"), "Metabolomics: Binary Proportion Distribution")
save_violin_plot(metabolomics_continuous, file.path(plot_dir, "Metabolomics_Continuous_ViolinPlot.pdf"), "Metabolomics: Continuous Proportion Distribution")

save_violin_plot(taxa_binary, file.path(plot_dir, "Taxa_Binary_ViolinPlot.pdf"), "Taxa: Binary Proportion Distribution")
save_violin_plot(taxa_continuous, file.path(plot_dir, "Taxa_Continuous_ViolinPlot.pdf"), "Taxa: Continuous Proportion Distribution")

save_violin_plot(concatenated_binary, file.path(plot_dir, "Concatenated_Binary_ViolinPlot.pdf"), "Concatenated: Binary Proportion Distribution")
save_violin_plot(concatenated_continuous, file.path(plot_dir, "Concatenated_Continuous_ViolinPlot.pdf"), "Concatenated: Continuous Proportion Distribution")