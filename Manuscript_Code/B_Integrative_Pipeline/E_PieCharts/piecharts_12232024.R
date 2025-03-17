# Load necessary library
library(ggplot2)
library(dplyr)

setwd('/Users/suzettepalmer/Desktop/IntegratedLearner_Dec2024/Piecharts_Figure2')
# Define the function
export_pie_chart <- function(data, output_file = "pie_chart.pdf") {
  # Add the Proportion column
  data <- data %>%
    mutate(Proportion = Count / sum(Count)) %>%
    arrange(desc(Dataset)) %>%
    mutate(Midpoint = cumsum(Proportion) - Proportion / 2)
  
  # Create the pie chart
  pie_chart <- ggplot(data, aes(x = "", y = Proportion, fill = Dataset)) +
    geom_bar(stat = "identity", width = 1, color = "black") + # Add black border
    coord_polar("y", start = 0) +
    theme_void() +
    geom_text(
      aes(
        x = 1.3, # Position labels outside the pie
        y = Midpoint,
        label = paste0(Dataset, " (", round(Proportion * 100, 1), "%)")
      ),
      size = 4
    ) +
    labs(title = "Number of Models Ranked in Top 3 per Dataset") +
    scale_fill_brewer(palette = "Set3") # You can change the color palette here
  
  # Export to PDF
  ggsave(output_file, plot = pie_chart, width = 8, height = 6)
  
  # Notify user
  message("Pie chart saved to ", output_file)
}


#filtered_testing_rmse
filtered_testing_rmse_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(2, 3, 39, 6, 4)
)
export_pie_chart(filtered_testing_rmse_integration, "filtered_testing_rmse_integration.pdf")

filtered_testing_rmse_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(13, 24, 17)
)
export_pie_chart(filtered_testing_rmse_model, "filtered_testing_rmse_model.pdf")

#filtered_testing_aucroc
filtered_testing_aucroc_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(6,7,7,8,5)
)
export_pie_chart(filtered_testing_aucroc_integration, "filtered_testing_aucroc_integration.pdf")

filtered_testing_aucroc_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(0,28,5)
)
export_pie_chart(filtered_testing_aucroc_model, "filtered_testing_aucroc_model.pdf")

###########
#Unfiltered_testing_rmse
unfiltered_testing_rmse_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(0,4,41,5,4)
)
export_pie_chart(unfiltered_testing_rmse_integration, "unfiltered_testing_rmse_integration.pdf")

unfiltered_testing_rmse_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(13, 25,16)
)
export_pie_chart(unfiltered_testing_rmse_model, "unfiltered_testing_rmse_model.pdf")

#filtered_testing_aucroc
unfiltered_testing_aucroc_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(5,10, 6,7,5)
)
export_pie_chart(unfiltered_testing_aucroc_integration, "unfiltered_testing_aucroc_integration.pdf")

unfiltered_testing_aucroc_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(2, 24, 7)
)
export_pie_chart(unfiltered_testing_aucroc_model, "unfiltered_testing_aucroc_model.pdf")

###########
#filtered_training_rmse
filtered_training_rmse_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(0, 0, 22, 5, 27)
)
export_pie_chart(filtered_training_rmse_integration, "filtered_training_rmse_integration.pdf")

filtered_training_rmse_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(14, 23, 17)
)
export_pie_chart(filtered_training_rmse_model, "filtered_training_rmse_model.pdf")

#filtered_training_aucroc
filtered_training_aucroc_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(0, 4, 9, 10, 13)
)
export_pie_chart(filtered_training_aucroc_integration, "filtered_training_aucroc_integration.pdf")

filtered_training_aucroc_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(7, 20, 9)
)
export_pie_chart(filtered_training_aucroc_model, "filtered_training_aucroc_model.pdf")

###########
#unfiltered_training_rmse
unfiltered_training_rmse_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(0,0,21,6,27)
)
export_pie_chart(unfiltered_training_rmse_integration, "unfiltered_training_rmse_integration.pdf")

unfiltered_training_rmse_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(14, 22, 18)
)
export_pie_chart(unfiltered_training_rmse_model, "unfiltered_training_rmse_model.pdf")

#unfiltered_testing_aucroc
unfiltered_training_aucroc_integration <- data.frame(
  Dataset = c("Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS"),
  Count = c(1, 3, 9, 10,13)
)
export_pie_chart(unfiltered_training_aucroc_integration, "unfiltered_training_aucroc_integration.pdf")

unfiltered_training_aucroc_model <- data.frame(
  Dataset = c("Elastic Net", "Random Forest", "XGBoost"),
  Count = c(14,7,15)
)
export_pie_chart(unfiltered_training_aucroc_model, "unfiltered_training_aucroc_model.pdf")


