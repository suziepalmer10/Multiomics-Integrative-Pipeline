#SECTION 1: Load in necessary libraries. 
#note: libraries specific for functions will be stored in the function files. 
library(caret)
library(kernlab)
library(xgboost)
library(tidyverse)
library(MLmetrics)
library(nnls)
library(optparse)
library(argparse)
library(glmnet)
library(pls)
library(pROC)
library(irr)
library(glmnet)
library(randomForest)
library(xgboost)

#set seed
set.seed(123)

#SECTION 1: DEBUGGING SECTION (OPTIONAL)
# Seed, fold/repeat and training/testing proportions
# Create a parser object
parser <- ArgumentParser(description = "Process integrated learner arguments")
# Define the arguments
parser$add_argument("--model_to_run", type = "character", help = "Model to Run")
parser$add_argument("--file_path", type = "character", help = "Pathway to Integrated Learner Directory.")
parser$add_argument("--input_file", type = "character", help = "Name of Input File.")
# studyauthor_typeofrun_responsevariable_stratifyvariable_modelran_date
parser$add_argument("--study_name", type = "character", help = "Name of Study for Output Files.")
parser$add_argument("--type_of_analysis", type = "character", help = "Type of Analysis: binary or continuous.")
parser$add_argument("--response_variable", type = "character", help = "Name of Response Variable Column.")
parser$add_argument("--stratify_variable", type = "character", default = '', help = "Name of Stratify Variable Column. Default is ''.")
parser$add_argument("--training_proportion", type = "numeric", default = 0.8, help = "Training Proportion to use. Default is 0.8.")
parser$add_argument("--num_repeats", type = "numeric", default=3, help = "Number of repeats to use for Repeated v-fold CV. Default is 3.")
parser$add_argument("--num_folds", type = "numeric", default = 5, help = "Number of folds to use for Repeated v-fold CV. Default is 5.")
# Parse the arguments
args <- parser$parse_args()

# Access the arguments
file_path <- args$file_path
setwd(file_path)
model_to_run <- args$model_to_run
source(model_to_run)
input_file <- args$input_file
study_name <- args$study_name
type_of_analysis <- args$type_of_analysis
response_variable <- args$response_variable
stratify_variable <- args$stratify_variable
training_proportion <- args$training_proportion
num_repeats <- args$num_repeats
num_folds <- args$num_folds

while (sink.number() > 0) {
  sink()
}

# Start capturing the output to the file
output_file <- paste0("results/AnalysisRunDocumentation/", study_name, ".txt")
# Check if the file already exists
if (file.exists(output_file)) {
  # Remove the existing file
  file.remove(output_file)
  print(paste("Existing file", output_file, "has been removed."))
}
sink(output_file)

#file metadata check 
print(paste('The model ran for this analysis is: ', model_to_run))
print(paste('The filepath for this analysis is: ', file_path))
print(paste('The input file for this analysis is: ', input_file))
print(paste('The study name for this analysis is: ', study_name))
print(paste('The type of performance metrics for this analysis is: ', type_of_analysis))
print(paste('The name of the response variable for this anlaysis is: ', response_variable))
print(paste('The name of the stratify variable, if used, for this analysis is: ', stratify_variable))
print(paste('The proportion of data used for training is: ', training_proportion))
print(paste('The number of repeats used for Repeated V-fold cross validation is: ', num_repeats))
print(paste('The number of folds used for Repeated V-fold cross validation is: ', num_folds))

# Capture the start time
start_time <- Sys.time()
# Define the metrics for continuous and binary data

#SECTION 2: LOAD IN DATA AND CHECK DATA TYPE AND LENGTH
data = read_csv(input_file, show_col_types=FALSE)
#checks
print(paste('The data file used is: ', input_file))
print(paste('The proportion used for testing is: ', training_proportion))
print(paste('The number of repeats is: ', num_repeats))
print(paste('The number of folds is: ', num_folds))
#DATA CHECK FOR METABOLITES
# Get the column names that start with 'm__'
m_columns <- grep("^m__", colnames(data), value = TRUE)
print(paste('The number of metabolites are: ', length(m_columns)))
# Check for missing values in metabolite columns
missing_values_m <- sapply(data[m_columns], function(column) any(is.na(column)))
columns_with_na_m <- length(names(missing_values_m)[missing_values_m])
print(paste("Metabolite number of columns with NA: ", columns_with_na_m))
#DATA CHECK FOR TAXA
#Get the column names that start with 't__'
t_columns <- grep("^t__", colnames(data), value = TRUE)
print(paste('The number of taxa are: ', length(t_columns)))
# Check for missing values in taxa columns
missing_values_t <- sapply(data[t_columns], function(column) any(is.na(column)))
columns_with_na_t <- length(names(missing_values_t)[missing_values_t])
print(paste("Taxa number of columns with NA: ", columns_with_na_t))
#DATA CHECK FOR RESPONSE VARIABLE AND STRATIFIED VARIABLE 
print(paste('TRUE/FALSE: The response variable, ', response_variable, ' has missing values. ', any(is.na(data[[response_variable]]))))
#check if stratify variable is present
# Check if the input is empty and assign NULL if so
if (stratify_variable == "") {
  stratify_variable <- NULL
  print('Notice: The stratify variable is NULL. ')
}
#if the stratify_variable is not NULL, check for missing values. 
if (!is.null(stratify_variable)) {
  print(paste('TRUE/FALSE: The stratify variable, ', stratify_variable, ' has missing values. ', any(is.na(data[[stratify_variable]]))))
}

# Check for missing values in the specified columns
#DATA CHECK BEFORE READING IN ANYMORE FILES - CHECK FOR NA
check_missing_values <- function(column_names) {
  any(sapply(data[column_names], function(column) any(is.na(column))))
}
# Combine all columns to check
columns_to_check <- c(response_variable, m_columns, t_columns)
#if stratify variable is not null
if (!is.null(stratify_variable)) {
  columns_to_check <- c(columns_to_check, stratify_variable)
}
# Check for missing values
if (check_missing_values(columns_to_check)) {
  stop("There are NA values in the dataset. Please check the data again before proceeding.")
} else {
  print("No missing values detected in the dataset. Analysis will proceed.")
}

# Check and convert columns to numeric if they are not already numeric
convert_to_numeric <- function(column_names) {
  for (col in column_names) {
    if (!is.numeric(data[[col]])) {
      # Attempt to convert to numeric
      data[[col]] <- as.numeric(data[[col]])
      if (any(is.na(data[[col]]))) {
        warning(paste("Some values in column", col, "could not be converted to numeric and have been set to NA."))
      }
    }
  }
}

# Convert specified columns to numeric
convert_to_numeric(c(m_columns, t_columns))
# Check if conversion was successful
check_numeric <- function(column_names) {
  all(sapply(data[column_names], is.numeric))
}
if (!check_numeric(c(m_columns, t_columns))) {
  stop("Not all metabolite and taxa columns are numeric. Please check the data.")
} else {
  print("All metabolite and taxa columns are numeric.")
}
# Handle script reading and response variable conversion
if (type_of_analysis == "binary") {
  #load in binary_functions for downstream analysis
  #these include accuracy, aucroc, and kappa.
  source("model_functions/binary_functions.R")
  binary_metrics <- list(
    accuracy = accuracy_calculation,
    kappa = kappa_calculation,
    aucroc = aucroc_calculation
  )
  metrics <- binary_metrics
  print("For binary table this is the counts recorded: ")
  print(table(data[[response_variable]]))
  print('Binary functions Accuracy, AUCROC and Kappa have been loaded.')
  # Convert response variable to factor for binary analysis
  if (!is.factor(data[[response_variable]])) {
    data[[response_variable]] <- as.factor(data[[response_variable]])
    print(paste("The response variable", response_variable, "has been converted to a factor."))
  }
} else if (type_of_analysis == "continuous") {
  #load in binary_functions for downstream analysis
  #these include mae, rmse and r^2
  source("model_functions/continuous_functions.R")
  continuous_metrics <- list(
    rmse = rmse_calculation,
    r2 = r2_calculation,
    mae = mae_calculation
  )
  metrics <- continuous_metrics
  print('Continuous functions R^2, RMSE and MAE have been loaded.')
  # Ensure response variable is numeric for continuous analysis
  if (!is.numeric(data[[response_variable]])) {
    data[[response_variable]] <- as.numeric(data[[response_variable]])
    print(paste("The response variable", response_variable, "has been converted to numeric."))
  }
} else {
  stop("Invalid Type of analysis. Please enter 'binary' or 'continuous'.")
}

#SECTION 3 - TRAINING AND TESTING BASED ON CONTINUOUS OR BINARY FUNCTIONS.
# Split train and test data
if (!is.null(stratify_variable)) {
  parts <- createDataPartition(data[[stratify_variable]], p=training_proportion, list=FALSE)
} else {
  parts <- createDataPartition(data, p=training_proportion, list=FALSE)
}
train <- data[parts, ]
test <- data[-parts, ]

# Storage for validation predictions
stacked_metabolomics_predictions <- list()
stacked_mss_predictions <- list()
ground_truth_predictions <- list()

#storage of testing predictions 
test_stacked_metabolomics_predictions <- list()
test_stacked_mss_predictions <- list()
test_stacked_response_predictions <- list()

#Storage for models created. 
metab_model <- list()
mss_model <- list()
concat_model <- list()

# #Subset the for concatenation, metabolomics and mss
# testset_concatenation <- subset(test, select = c(response_variable, m_columns, t_columns))
# testset_metabolomics <- subset(test, select = c(response_variable, m_columns))
# testset_mss <- subset(test, select = c(response_variable, t_columns))

# Subset the training data
testset_concatenation  <- test %>% select(response_variable, m_columns, t_columns)
testset_metabolomics <- test %>% select(response_variable, m_columns)
testset_mss<- test %>% select(response_variable, t_columns)


# Perform repeated cross-validation
for (repeat_ in 1:num_repeats) {
  for (fold in 1:num_folds) {
    # Create train and test indices for the current fold
    set.seed(repeat_ * fold)  # Vary the seed for each fold
    indices <- createFolds(train[[response_variable]], k = num_folds, list = TRUE)
    train_index <- unlist(indices[-fold])
    validation_index <- unlist(indices[fold])
    
    # Split the data into training and testing sets
    train_data <- train[train_index, ]
    validation_data <- train[validation_index, ]
    
    # Subset the training data
    train_concat <- train_data %>% select(response_variable, m_columns, t_columns)
    train_metabolomics <- train_data %>% select(response_variable, m_columns)
    train_mss <- train_data %>% select(response_variable, t_columns)
    # Subset the validation data
    validation_concat <- validation_data %>% select(response_variable, m_columns, t_columns)
    validation_metabolomics <- validation_data %>% select(response_variable, m_columns)
    validation_mss <- validation_data %>% select(response_variable, t_columns)
    
    # Models to run
    metab_predictions <- model(train_metabolomics, validation_metabolomics, testset_metabolomics, response_variable, type_of_analysis)
    mss_predictions <- model(train_mss, validation_mss, testset_mss, response_variable, type_of_analysis)
    concat_predictions <- model(train_concat, validation_concat, testset_concatenation, response_variable, type_of_analysis)
    
    # Store models to revisit
    metab_model[[length(metab_model)+1]] <- metab_predictions[[5]]
    mss_model[[length(mss_model)+1]] <- mss_predictions[[5]]
    concat_model[[length(concat_model)+1]] <- concat_predictions[[5]]
    
    # Store predictions for downstream stacking analysis
    # This is for the validation set
    stacked_metabolomics_predictions[[length(stacked_metabolomics_predictions)+1]] <- metab_predictions[[1]]
    stacked_mss_predictions[[length(stacked_mss_predictions)+1]] <- mss_predictions[[1]]
    ground_truth_predictions[[length(ground_truth_predictions)+1]] <- metab_predictions[[2]]
    # This is for the testing set
    test_stacked_metabolomics_predictions[[length(test_stacked_metabolomics_predictions)+1]] <- metab_predictions[[3]]
    test_stacked_mss_predictions[[length(test_stacked_mss_predictions)+1]] <- mss_predictions[[3]]
    test_stacked_response_predictions[[length(test_stacked_response_predictions)+1]] <- metab_predictions[[4]]
    
    # Calculate and store the fold-level performance metrics
    for (metric_name in names(metrics)) {
      metric_function <- metrics[[metric_name]]
      assign(paste0(metric_name, "_fold_metabolomics"), c(get(paste0(metric_name, "_fold_metabolomics")), metric_function(metab_predictions[[2]], metab_predictions[[1]])))
      assign(paste0(metric_name, "_fold_mss"), c(get(paste0(metric_name, "_fold_mss")), metric_function(mss_predictions[[2]], mss_predictions[[1]])))
      assign(paste0(metric_name, "_fold_concat"), c(get(paste0(metric_name, "_fold_concat")), metric_function(concat_predictions[[2]], concat_predictions[[1]])))
      assign(paste0("test_", metric_name, "_fold_metabolomics"), c(get(paste0("test_", metric_name, "_fold_metabolomics")), metric_function(metab_predictions[[4]], metab_predictions[[3]])))
      assign(paste0("test_", metric_name, "_fold_mss"), c(get(paste0("test_", metric_name, "_fold_mss")), metric_function(mss_predictions[[4]], mss_predictions[[3]])))
      assign(paste0("test_", metric_name, "_fold_concat"), c(get(paste0("test_", metric_name, "_fold_concat")), metric_function(concat_predictions[[4]], concat_predictions[[3]])))
    }
  }
}

#SECTION 4: Integrated performance models 

#TRAINING DATA
#loop through the stacked metabolomics, since this list will be the same size 
#as mss and ground truth. 

# Lasso weights and lambda
sparse_nnls_weights <- list()
best_lambda_list <- list()
best_alpha_list <- list()
# NNLS weights
nnls_weights <- list()
# PLS models
pls_models <- list()
best_ncomp_list <- list()

for (i in 1:length(stacked_metabolomics_predictions)) {
  New_df_train <- cbind(
    unlist(stacked_metabolomics_predictions[[i]]), 
    unlist(stacked_mss_predictions[[i]]), 
    unlist(ground_truth_predictions[[i]])
  )
  colnames(New_df_train) <- c('Metabolomics', 'MSS', 'Ground Truth')
  New_df_train <- as.data.frame(New_df_train)
  
  # NNLS weighted
  train_weights <- nnls(as.matrix(New_df_train[, 1:2]), as.matrix(New_df_train[, 3]))
  New_df_train$weighted <- New_df_train$Metabolomics * train_weights$x[1] + New_df_train$MSS * train_weights$x[2]
  # Averaged stacked
  New_df_train$average <- New_df_train$Metabolomics * 0.5 + New_df_train$MSS * 0.5
  
  # Determine family for glmnet
  if (type_of_analysis == "binary") {
    family <- "binomial"
  } else if (type_of_analysis == "continuous") {
    family <- "gaussian"
  } else {
    stop("Unsupported response type. Please use 'binary' or 'continuous'.")
  }
  
  # Lasso with non-negative lower bounds
  alpha_values <- seq(0, 1, by = 0.1)
  lambda_values <- 10^seq(2, -2, by = -0.1)
  best_alpha <- 0
  best_lambda <- 0
  best_cv_error <- Inf
  
  for (alpha in alpha_values) {
    cv_fit <- cv.glmnet(
      as.matrix(New_df_train[, 1:2]), 
      New_df_train$`Ground Truth`, 
      alpha = alpha, 
      nfolds = 5, 
      family = family,
      lambda = lambda_values
    )
    cv_error <- min(cv_fit$cvm)
    
    if (cv_error < best_cv_error) {
      best_cv_error <- cv_error
      best_alpha <- alpha
      best_lambda <- cv_fit$lambda.min
    }
  }
  
  final_fit <- glmnet(
    as.matrix(New_df_train[, 1:2]), 
    New_df_train$`Ground Truth`, 
    alpha = best_alpha, 
    lambda = best_lambda, 
    family = family
  )
  sparse_weighted_predictions <- as.data.frame(predict(final_fit, as.matrix(New_df_train[, 1:2])))
  New_df_train$sparse_weighted <- sparse_weighted_predictions$s0
  
  # PLS model fitting and cross-validation
  pls_cv <- plsr(`Ground Truth` ~ Metabolomics + MSS, data = New_df_train, validation = "CV")
  best_ncomp <- which.min(pls_cv$validation$PRESS)
  pls_model <- plsr(`Ground Truth` ~ Metabolomics + MSS, data = New_df_train, ncomp = best_ncomp)
  pls_models[[i]] <- pls_model
  best_ncomp_list[[i]] <- best_ncomp
  
  #New_df_train$pls_pred <- predict(pls_model, New_df_train[, 1:2], ncomp = best_ncomp)
  pls_train_prediction <- as.data.frame(predict(pls_model, New_df_train[, c("Metabolomics", "MSS")], ncomp = best_ncomp))
  New_df_train$pls_pred <- drop(pls_train_prediction[, 1])
  
  for (metric_name in names(metrics)) {
    metric_function <- metrics[[metric_name]]
    assign(paste0(metric_name, "_averaged_stacked"), c(get(paste0(metric_name, "_averaged_stacked")), metric_function(New_df_train$`Ground Truth`, New_df_train$average)))
    assign(paste0(metric_name, "_sparse_nnls"), c(get(paste0(metric_name, "_sparse_nnls")), metric_function(New_df_train$`Ground Truth`, New_df_train$weighted)))
    assign(paste0(metric_name, "_weighted_nnls"), c(get(paste0(metric_name, "_weighted_nnls")), metric_function(New_df_train$`Ground Truth`, New_df_train$sparse_weighted)))
    assign(paste0(metric_name, "_pls"), c(get(paste0(metric_name, "_pls")), metric_function(New_df_train$`Ground Truth`, New_df_train$pls_pred)))
  }
  
  # Store NNLS weights and PLS model for testing
  nnls_weights[[i]] <- train_weights$x
  sparse_nnls_weights[[i]] <- coef(final_fit)
  best_lambda_list[[i]] <- best_lambda
  best_alpha_list[[i]] <- best_alpha
}

# TESTING DATA
test_sparse_nnls_weights <- list()
test_pls_models <- list()
new_df_test_list <- list()

for (i in 1:length(test_stacked_metabolomics_predictions)) {
  new_df_test <- cbind(
    unlist(test_stacked_metabolomics_predictions[[i]]), 
    unlist(test_stacked_mss_predictions[[i]]), 
    unlist(test_stacked_response_predictions[[i]])
  )
  colnames(new_df_test) <- c('Metabolomics', 'MSS', 'Ground Truth')
  new_df_test <- as.data.frame(new_df_test)
  new_df_test$weighted <- new_df_test$Metabolomics * nnls_weights[[i]][1] + new_df_test$MSS * nnls_weights[[i]][2]
  new_df_test$average <- new_df_test$Metabolomics * 0.5 + new_df_test$MSS * 0.5
  
  # Lasso/enet using glmnet with the best lambda from training
  final_fit <- glmnet(as.matrix(new_df_test[, 1:2]), new_df_test$`Ground Truth`, family = family, alpha = best_alpha_list[[i]], lambda = best_lambda_list[[i]])
  sparse_test_weighted_predictions <- as.data.frame(predict(final_fit, as.matrix(new_df_test[, 1:2])))
  new_df_test$sparse_weighted <- sparse_test_weighted_predictions$s0
  test_sparse_nnls_weights[[i]] <- coef(final_fit)
  
  # Predict using PLS model
  pls_model <- pls_models[[i]]
  #new_df_test$pls_pred <- predict(pls_model, new_df_test[, 1:2], ncomp = best_ncomp_list[[i]])
  pls_test_predictions <- as.data.frame(predict(pls_model, new_df_test[, c("Metabolomics", "MSS")], ncomp = best_ncomp_list[[i]]))
  new_df_test$pls_pred <- unlist(pls_test_predictions[, 1])
  new_df_test_list[[i]] <-  new_df_test
    
  for (metric_name in names(metrics)) {
    metric_function <- metrics[[metric_name]]
    assign(paste0('test_', metric_name, "_averaged_stacked"), c(get(paste0('test_', metric_name, "_averaged_stacked")), metric_function(new_df_test$`Ground Truth`, new_df_test$average)))
    assign(paste0('test_', metric_name, "_sparse_nnls"), c(get(paste0('test_', metric_name, "_sparse_nnls")), metric_function(new_df_test$`Ground Truth`, new_df_test$weighted)))
    assign(paste0('test_', metric_name, "_weighted_nnls"), c(get(paste0('test_', metric_name, "_weighted_nnls")), metric_function(new_df_test$`Ground Truth`, new_df_test$sparse_weighted)))
    assign(paste0('test_', metric_name, "_pls"), c(get(paste0('test_', metric_name, "_pls")), metric_function(new_df_test$`Ground Truth`, new_df_test$pls_pred)))
  }
}

testing <- new_df_test_list[[2]]

lapply(testing, class)

#SECTION 5 - Results Analysis and Data Visualization
# Function to calculate means and standard deviations and return a formatted dataframe
stats_results <- function(lists, row_labels, metric) {
  means <- sapply(lists, function(x) mean(unlist(x), na.rm = TRUE))
  sds <- sapply(lists, function(x) sd(unlist(x), na.rm = TRUE))
  result_df <- data.frame(
    Category = row_labels,
    Mean = means,
    SD = sds
  )
  colnames(result_df)[2:3] <- paste(metric, c("_Mean", "_SD"), sep = " ")
  return(result_df)
}

# load in source file functions for correct visualizations
if (type_of_analysis == "binary") {
  source('model_functions/binary_performance_functions.R')
} else if (type_of_analysis == "continuous") {
  source('model_functions/continuous_performance_functions.R')
} else {
  stop("Unsupported response type. Please use 'binary' or 'continuous'.")
}

# Convert the data to long format
long_data <- all_results_df %>%
  pivot_longer(cols = -Category, names_to = "Metric", values_to = "Value") %>%
  # Split column names into Metric and Type based on the format "Metric Type"
  separate(Metric, into = c("Metric", "Type"), sep = "_") %>%
  # Ensure Type is treated as a factor with levels "Mean" and "SD"
  mutate(Type = factor(Type, levels = c("Mean", "SD")))

# Separate Mean and SD for plotting
mean_df <- long_data %>% filter(Type == "Mean")
sd_df <- long_data %>% filter(Type == "SD")

plot_df <- mean_df %>%
  left_join(sd_df, by = c("Category", "Metric"), suffix = c(".mean", ".sd")) %>%
  mutate(
    ymin = Value.mean - Value.sd,
    ymax = Value.mean + Value.sd
  )

# Update plot_df with ordered factors
plot_df_order <- plot_df %>%
  mutate(
    Metric = factor(Metric, levels = facet_order),
    Category = factor(Category, levels = x_order)
  )

#create a faceted boxplot for the performance metrics 
plot <- ggplot(data = plot_df_order, aes(x = Category, y = Value.mean, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  labs(title = study_name, x = "Integration Technique", y = "Metric Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#SECTION 6: SAVE RESULTS 
#save the plot file 
plot_file = paste0("results/visualizations/", study_name, ".png")
ggsave(filename = plot_file, plot = plot, width = 10, height = 8, dpi = 300, bg='white')
#save the RData file
rdata_filename = paste0("results/RDataFiles/", study_name, ".RData")
save.image(file=rdata_filename)
#save the performance metrics file
performance_filename = paste0("results/performance_csv/", study_name, ".csv")
write.csv(all_results_df, file = performance_filename, row.names = TRUE)


# Capture the end time
end_time <- Sys.time()
# Calculate the elapsed time
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
# Convert elapsed time to hours and minutes
hours <- floor(elapsed_time / 3600)
minutes <- floor((elapsed_time %% 3600) / 60)
# Print the elapsed time
cat(sprintf("Elapsed time: %d hours and %d minutes or %f seconds\n", hours, minutes, elapsed_time))
sink()
#############

