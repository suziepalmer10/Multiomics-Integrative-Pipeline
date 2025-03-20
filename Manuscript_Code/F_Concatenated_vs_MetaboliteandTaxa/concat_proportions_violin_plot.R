library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

#can be full_dimensional or feature_reduced
#run_info = "feature_reduced"
run_info = "full_dimensional"

# Define directories and corresponding file names
top_levels <- c(1, 5, 10, 20, 50)
base_dir <- "/Users/suzettepalmer/Desktop/Integrated_Pipeline_March2025/Code_and_Analyses/Manuscript_Code/F_Concatenated_vs_MetaboliteandTaxa/"  # Set your base directory

# Define the sheet names
sheets <- c("enet", "rf", "xgb")

# Initialize a list to store the data
data_list <- list()

# Loop through each directory
for (n in top_levels) {
  dir_path <- file.path(base_dir, paste0("Top", n))
  if (run_info == "feature_reduced") {
    file_name <- paste0(n, "featurereduced_concat_proportions.xlsx")
  }
  if (run_info == "full_dimensional") {
    file_name <- paste0(n, "fulldimensional_concat_proportions.xlsx")
  } 
  
  file_path <- file.path(dir_path, file_name)
  
  # Check if file exists before reading
  if (file.exists(file_path)) {
    # Read each sheet into a list
    sheet_data <- lapply(sheets, function(sheet) {
      read_excel(file_path, sheet = sheet)
    })
    names(sheet_data) <- sheets
    
    # Store in the main list
    data_list[[paste0("Top", n)]] <- sheet_data
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Initialize an empty dataframe to store combined data
plot_data <- data.frame()

# Loop through each top level dataset
for (top_level in names(data_list)) {
  # Extract the numeric value from the "Top" label (e.g., "Top5" -> 5)
  top_n <- as.numeric(gsub("Top", "", top_level))
  
  for (model in sheets) {
    if (!is.null(data_list[[top_level]][[model]])) {
      df <- data_list[[top_level]][[model]]
      
      # Check if required columns exist
      if ("Taxa" %in% colnames(df) & "Metabolites" %in% colnames(df) & "row_id" %in% colnames(df)) {
        
        # Compute proportions while keeping row_id intact
        df_long <- df %>%
          mutate(across(c(Taxa, Metabolites), ~ . / top_n)) %>%  # Divide only Taxa and Metabolites
          pivot_longer(cols = c(Taxa, Metabolites), names_to = "Feature_Type", values_to = "Proportion") %>%
          mutate(Model = model, Top_Level = top_level) %>%
          select(row_id, Feature_Type, Proportion, Model, Top_Level)  # Ensure row_id is retained
        
        # Append to main data
        plot_data <- bind_rows(plot_data, df_long)
      }
    }
  }
}


fd_continuous <- c("nf_era_age" , "nf_era_cholesterol", "nf_era_glucose",
                   "un_franzosa_ControlvsCD_Age_11122024", "un_franzosa_ControlvsCD_Fp_11122024",
                   "un_franzosa_ControlvsDisease_Age_11122024", "un_franzosa_ControlvsDisease_Fp_11122024",
                   "un_franzosa_ControlvsUC_Age_11122024", "un_franzosa_ControlvsUC_Fp_11122024",
                   "nf_wang_age_12092024", "nf_wang_bmi_12092024", "nf_wang_creatinine_12092024",
                   "nf_wang_egfr_12092024", "nf_wang_urea_12092024", "nf_yachida_age_12092024",
                   "nf_yachida_alcohol_12092024", "nf_yachida_BrinkmanIndex_12092024",
                   "nf_yachida_BMI_12092024" )
fd_binary <- c("nf_era_alcohol",  "nf_era_gender", "nf_era_sg",
               "un_franzosa_ControlvsCD_ConvCD_11122024", "un_franzosa_ControlvsDisease_ConvDisease_11122024",
               "un_franzosa_ControlvsUC_ConvUC_11122024", "nf_wang_studygroup_12092024",
               "nf_yachida_gender_12092024", "nf_yachida_healthyvscancer_12092024",
               "nf_yachida_healthyvsstageIII_IV_12092024" , "nf_wang_gender_12092024",
               "nf_yachida_healthyvsearly_12092024", "nf_yachida_healthyvsstageI_II_12092024")



binary_filter <- c('era_alcohol', "era_sg", "franzosa_ControlvsCD_ConvCD_11122024",
                   "franzosa_ControlvsDisease_ConvDisease_11122024", "franzosa_ControlvsUC_ConvUC_11122024",
                   "wang_studygroup_12092024", "yachida_gender_12092024", "yachida_healthyvscancer_12092024",
                   "yachida_healthyvsearly_12092024", "yachida_healthyvsstageI_II_12092024",
                   "yachida_healthyvsstageIII_IV_12092024", "era_gender", "wang_gender_12092024" )
continuous_filter <- c('era_age', 'era_cholesterol', "era_glucose", "franzosa_ControlvsCD_Age_11122024",
                       "franzosa_ControlvsCD_Fp_11122024", "franzosa_ControlvsDisease_Age_11122024",
                       "franzosa_ControlvsDisease_Fp_11122024", "franzosa_ControlvsUC_Age_11122024",
                       "franzosa_ControlvsUC_Fp_11122024", "wang_age_12092024", "wang_creatinine_12092024",
                       "wang_egfr_12092024", "wang_urea_12092024", "yachida_age_12092024", "yachida_alcohol_12092024",
                       "yachida_BrinkmanIndex_12092024", "wang_bmi_12092024", "yachida_BMI_12092024")

plot_data$row_id <- gsub("(_enet_concat_importance|_xgb_concat_importance|_xgboost_concat_importance|_rf_concat_importance)$", "", plot_data$row_id)

custom_colors <- c("Metabolites" = '#780C28',
            "Taxa" = "#0000FF")



violin_plot <- function(df, filter_select, title_) { 
df$Top_Level <- factor(df$Top_Level, levels = c("Top1", "Top5", "Top10", "Top20", "Top50"))
df_filtered <- df %>% filter(row_id %in% filter_select)
# Create violin plot with proportions
g<- ggplot(df_filtered, aes(x = Feature_Type, y = Proportion, fill = Feature_Type)) +
  geom_violin(alpha = 0.7, trim = TRUE, width = 0.6) +
  geom_jitter(aes(color = Model), position = position_jitter(width = 0.2, height = 0), 
              alpha = 0.6, size = 1.5, shape = 1, color = 'black') + 
  facet_grid(Top_Level ~ Model) +
  labs(title = title_,
       x = "Feature Type",
       y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors for boxplot
  #scale_fill_manual(values = c("Taxa" = "blue", "Metabolites" = "red"))
  scale_y_continuous(limits = c(0, 1.0), expand = expansion(mult = c(0.05, 0.05))) 
} 


generate_violin_plot <- function(df, filter, title_cont) {
  violin <- violin_plot(df, filter, title_cont)
  pdf(paste0(base_dir, title_cont, "_violin_plot.pdf"), width = 6, height = 8)
  print(violin)
  dev.off()
}

if (run_info == 'feature_reduced') {
  generate_violin_plot(plot_data, continuous_filter,"FeatureRed_Concat_FeatureProportions_Continuous")
  generate_violin_plot(plot_data, binary_filter,"FeatureRed_Concat_FeatureProportions_Binary")
}

if (run_info == "full_dimensional") {
  generate_violin_plot(plot_data, fd_continuous,"FullDimensional_Concat_FeatureProportions_Continuous")
  generate_violin_plot(plot_data, fd_binary,"FullDimensional_Concat_FeatureProportions_Binary")
}


