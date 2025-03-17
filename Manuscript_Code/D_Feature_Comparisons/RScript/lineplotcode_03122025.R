# Load the library
library(readxl)
library(ggplot2)
library(tidyverse)

setwd('/Users/suzettepalmer/Desktop/IntegratedLearner_March2025/Code_and_Files/Feature_Comparisons')

top1_df <- read_excel('Top1/1metab_filtered_03112025.xlsx')
top3_df <- read_excel('Top3/3metab_filtered_03112025.xlsx')
top5_df <- read_excel('Top5/5metab_filtered_03112025.xlsx')
top10_df <- read_excel('Top10/10metab_filtered_03112025.xlsx')
top20_df <- read_excel('Top20/20metab_filtered_03112025.xlsx')
top50_df <- read_excel('Top50/50metab_filtered_03112025.xlsx')


# Extract the DataFrame_Names reference from top1_df
valid_names <- top1_df$Dataframe_Name

# Define a function to align dataframe names
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

# Store all datasets in a named list
data_list <- list(
  "Top1" = top1_df,
  "Top3" = top3_df,
  "Top5" = top5_df,
  "Top10" = top10_df,
  "Top20" = top20_df,
  "Top50" = top50_df
)

# Process data into long format while ensuring missing dataframe handling
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


custom_colors <- c("all_overlap" = "#1E88E5",  # Blue
                   "enet_rf_overlap" = "#D55E00",  # Orange
                   "rf_only" = "#117733",  # Green
                   "xgb_only" = "#666666",  # Red
                   "enet_only" = "#CCAA00", 
                   "enet_xgb_overlap" = "#332288", 
                   "rf_xgb_overlap" = "#604020")  # Purple
custom_shapes <- c("all_overlap" = 21,  # Square solid
                   "enet_rf_overlap" = 22,  # Circle open
                   "rf_only" = 21,  # Triangle down solid
                   "xgb_only" = 24,  # Square cross
                   "enet_only" = 25,   #Triangle up open 
                   "enet_xgb_overlap" = 23,   # Plus
                   "rf_xgb_overlap" = 21)   # Triangle Up solid
df_combined$Metric <- factor(df_combined$Metric, levels = names(custom_colors))
# Create the plot with facets (6 columns, 1 row)

pdf("my_plot.pdf", width = 8, height = 11)
ggplot(df_combined, aes(x = Dataframe_Name, y = Value, color = Metric, shape =  Metric, group = Metric)) +
  geom_line(size = .25, na.rm = TRUE) + 
  geom_point(size = 3, na.rm = TRUE) +
  theme_minimal() +
  labs(title = "Comparisons of Featured Reduced Metabolites Selected Across Datasets", 
       x = "Dataframe Name", y = "Proportion of Features") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),  # Keep y-axis text readable
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels vertically
    strip.text = element_text(size = 8),  # Formatting for facet titles
    panel.spacing = unit(0.05, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_wrap(~Dataset, nrow = 6, scales = "fixed")+
  scale_color_manual(values = custom_colors)+  # Apply custom colors
  scale_shape_manual(values = custom_shapes)+  # Apply custom shapes
  scale_y_continuous(limits = c(0, 1.2), expand = expansion(mult = c(0.2, 0.2))) 


dev.off()



