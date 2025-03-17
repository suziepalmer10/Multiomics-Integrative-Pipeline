# Load required libraries
library(readxl)
library(ggplot2)
library(tidyverse)

setwd('/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/03112025_Coefficients/')

top1_df <- read_excel('Top1/1metab_filtered_03112025.xlsx')
top3_df <- read_excel('Top3/3metab_filtered_03112025.xlsx')
top5_df <- read_excel('Top5/5metab_filtered_03112025.xlsx')
top10_df <- read_excel('Top10/10metab_filtered_03112025.xlsx')
top20_df <- read_excel('Top20/20metab_filtered_03112025.xlsx')
top50_df <- read_excel('Top50/50metab_filtered_03112025.xlsx')

# Extract the DataFrame_Names reference from top1_df
valid_names <- top1_df$Dataframe_Name

# Function to align dataframe names
align_dataframe_names <- function(df, valid_names) {
  missing_names <- setdiff(valid_names, df$Dataframe_Name)
  
  if (length(missing_names) > 0) {
    missing_rows <- data.frame(Dataframe_Name = missing_names, matrix(NA, nrow = length(missing_names), ncol = ncol(df) - 1))
    colnames(missing_rows) <- colnames(df)  # Ensure column names match
    df <- bind_rows(df, missing_rows)
  }
  df <- df %>% arrange(factor(Dataframe_Name, levels = valid_names))
  return(df)
}

# Apply function to all datasets
top3_df <- align_dataframe_names(top3_df, valid_names)
top5_df <- align_dataframe_names(top5_df, valid_names)
top10_df <- align_dataframe_names(top10_df, valid_names)
top20_df <- align_dataframe_names(top20_df, valid_names)
top50_df <- align_dataframe_names(top50_df, valid_names)

# Store datasets in a named list
data_list <- list(
  "Top1" = top1_df,
  "Top3" = top3_df,
  "Top5" = top5_df,
  "Top10" = top10_df,
  "Top20" = top20_df,
  "Top50" = top50_df
)

# Process data into long format
long_data_list <- list()
for (name in names(data_list)) {
  df <- data_list[[name]]
  if (!is.null(df)) {
    long_data_list[[name]] <- df %>%
      pivot_longer(cols = -Dataframe_Name, names_to = "Metric", values_to = "Value") %>%
      mutate(Dataset = name)
  } 
}

# Combine all datasets
df_combined <- bind_rows(long_data_list)

df_combined$Metric <- factor(df_combined$Metric, levels = c("all_overlap", "enet_rf_overlap", "enet_xgb_overlap",
                                                            "rf_xgb_overlap", "enet_only", "rf_only", "xgb_only"))
df_combined$Dataset <- factor(df_combined$Dataset, levels = c("Top1", "Top3", "Top5", "Top10", "Top20", "Top50"))

df_combined$Metric <- trimws(df_combined$Metric)

binary <- c("era_alcohol", 'era_sg',  "franzosa_ControlvsDisease_ConvDisease", 
            "yachida_gender", "yachida_healthyvscancer", "yachida_healthyvsearly",               
            "yachida_healthyvsstageI_II", "yachida_healthyvsstageIII_IV", "wang_studygroup" )
continuous <- c("era_age", "era_cholesterol", "era_glucose", "franzosa_ControlvsCD_Age",             
"franzosa_ControlvsCD_Fp", "franzosa_ControlvsDisease_Age", "franzosa_ControlvsDisease_Fp",         
"franzosa_ControlvsUC_Age", "franzosa_ControlvsUC_Fp", "yachida_age", "yachida_alcohol", 
"yachida_BrinkmanIndex", "wang_age", "wang_creatinine","wang_egfr", "wang_urea"        
                )
df_combined_1 <- df_combined %>% filter(Dataset != "Top3")
  
df_binary <- df_combined_1[df_combined_1$Dataframe_Name %in% binary, ]
# Subset continuous dataframe
df_continuous <- df_combined_1[df_combined_1$Dataframe_Name %in% continuous, ]




# Custom colors for boxplot
custom_colors <- c("all_overlap" = "#0071bc",  # Blue
                   "enet_rf_overlap" = "#939598",  # Orange
                   "rf_only" = "#B31E2D",  # Green
                   "xgb_only" = "#F3C366",  # Gray
                   "enet_only" = "#6b00b9", 
                   "enet_xgb_overlap" = "#76dcd2", 
                   "rf_xgb_overlap" = "#FFFFFF")  # Purple

range(df_binary$Value, na.rm = TRUE)

violin_plot <- function(title_select, df) {
# Save the plot as a PDF
#pdf("boxplot_plot.pdf", width = 2, height = 3)
#while (dev.cur() > 1) dev.off()
ggplot(df, aes(x = Metric, y = Value, fill = Metric)) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_violin(alpha = 0.7, trim = FALSE, width = 0.6) +
  geom_jitter(aes(color = Metric), position = position_jitter(width = 0.2, height = 0), 
               alpha = 0.6, size = 1.5, shape = 1, color = 'black') + 
  theme_minimal() +
  labs(title = title_select, 
       x = "Metric", y = "Proportion of Features") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.05, "lines"),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.y = element_blank(),
    #panel.grid.major.x = element_blank(),  # Removes vertical major grid lines
    #panel.grid.minor.x = element_blank()
  ) +
  facet_wrap(~Dataset, nrow = 6, scales = "fixed") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors for boxplot
  #scale_color_manual(values = custom_colors) +  # Apply same custom colors for jitter points
  scale_y_continuous(limits = c(0, 1.0), expand = expansion(mult = c(0.05, 0.05))) 
}

metab_filter_cont <- violin_plot("Comparisons of Feature Reduced Continuous Metabolites Selected Across Datasets",
                                 df_continuous)
metab_filter_binary <- violin_plot("Comparisons of Feature Reduced Binary Metabolites Selected Across Datasets",
                                   df_binary)

pdf("metabolite_filtered_continuous_violin_plot.pdf", width = 6, height = 8)
metab_filter_cont
dev.off()
pdf("metabolite_filtered_binary_violin_plot.pdf", width = 6, height = 8)
metab_filter_binary
dev.off()

