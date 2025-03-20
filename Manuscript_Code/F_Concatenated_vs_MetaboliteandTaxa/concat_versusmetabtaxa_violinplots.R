library(readxl)
library(ggplot2)
library(dplyr)



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
    file_name <- paste0(n, "feature_reduced_comparisons.xlsx")
  }
  if (run_info == "full_dimensional") {
    file_name <- paste0(n, "full_dimensional_comparisons.xlsx")
  } else {
    print('No file found or invalid input.')
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



# Initialize an empty dataframe to store results
metab_df <- data.frame()
taxa_df <- data.frame()
# Loop through each Top level in data_list
for (top in names(data_list)) {
  for (model in c("enet", "rf", "xgb")) {
    # Check if the column exists in the sheet before extracting
    if ("Metab_Concat Overlap" %in% colnames(data_list[[top]][[model]])) {
      temp_df <- data_list[[top]][[model]] %>%
        select("Prefix", "Metab_Concat Overlap") %>%
        mutate(TopN = top, Model = model)  # Add columns for Top level and model type
      # Combine results
      metab_df <- bind_rows(metab_df, temp_df)
    } 
    if ("Taxa_Concat Overlap" %in% colnames(data_list[[top]][[model]])) {
      temp_df <- data_list[[top]][[model]] %>%
        select("Prefix", "Taxa_Concat Overlap") %>%
        mutate(TopN = top, Model = model)  # Add columns for Top level and model type
      # Combine results
      taxa_df <- bind_rows(taxa_df, temp_df)
    }
  }
}


taxa_df$Prefix <- gsub("(_enet|_xgb|_xgboost|_rf)$", "", taxa_df$Prefix)
metab_df$Prefix <- gsub("(_enet|_xgb|_xgboost|_rf)$", "", metab_df$Prefix)
taxa_df <- taxa_df %>% rename(Overlap= `Taxa_Concat Overlap`)
metab_df <- metab_df %>% rename(Overlap= `Metab_Concat Overlap`)


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



violin_plot <- function(df, filter_select, title_) {

#filter binary or continuous features
df$TopN <- factor(df$TopN, levels = c("Top1", "Top5", "Top10", "Top20", "Top50"))
custom_colors <- c("enet" = "#FFFFFF", 
                   "rf" = "#939598", 
                   "xgb"= "#6b00b9")
df_filtered <- df %>% filter(Prefix %in% filter_select)

g <- ggplot(df_filtered, aes(x = Model, y = Overlap, fill = Model)) +
  #geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot with transparency
  #geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add jitter for visibility of points
  geom_violin(alpha = 0.7, trim = TRUE, width = 0.6) +
  geom_jitter(aes(color = Model), position = position_jitter(width = 0.2, height = 0), 
              alpha = 0.6, size = 1.5, shape = 1, color = 'black') + 
  theme_minimal() +
  labs(title = title_, 
       x = "Metric", y = "Proportion of Features") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.05, "lines"),
  ) +
  facet_wrap(~TopN, nrow = 6, scales = "fixed") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors for boxplot
  scale_y_continuous(limits = c(0, 1.0), expand = expansion(mult = c(0.05, 0.05))) 
  return(g)

}

generate_violin_plot <- function(df, filter, title_cont) {
  violin <- violin_plot(df, filter, title_cont)
  pdf(paste0(base_dir, title_cont, "_violin_plot.pdf"), width = 6, height = 8)
  print(violin)
  dev.off()
}

if (run_info == 'feature_reduced') {
  generate_violin_plot(metab_df, continuous_filter,"FeatureRed_Concat_Metab_Overlap_Continuous")
  generate_violin_plot(metab_df, binary_filter,"FeatureRed_Concat_Metab_Overlap_Binary")
  generate_violin_plot(taxa_df, continuous_filter,"FeatureRed_Concat_Taxa_Overlap_Continuous")
  generate_violin_plot(taxa_df, binary_filter,"FeatureRed_Concat_Taxa_Overlap_Binary")
}

if (run_info == "full_dimensional") {
  generate_violin_plot(metab_df, fd_continuous,"FullDimensional_Concat_Metab_Overlap_Continuous")
  generate_violin_plot(metab_df, fd_binary,"FullDimensional_Concat_Metab_Overlap_Binary")
  generate_violin_plot(taxa_df, fd_continuous,"FullDimensional_Concat_Taxa_Overlap_Continuous")
  generate_violin_plot(taxa_df, fd_binary,"FullDimensional_Concat_Taxa_Overlap_Binary")
}




 