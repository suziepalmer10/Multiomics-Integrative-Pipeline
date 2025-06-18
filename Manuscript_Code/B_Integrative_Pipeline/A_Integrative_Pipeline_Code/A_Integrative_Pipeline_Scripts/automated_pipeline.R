#SECTION 1: Load in necessary libraries. 
#note: libraries specific for functions will be stored in the function files. 
pkgs_needed <- c(
  "caret", "kernlab", "xgboost", "tidyverse", "MLmetrics",
  "nnls", "argparse", "glmnet", "pls", "pROC", "irr", "randomForest","readr"
)

# TRUE/FALSE test + action
for (pkg in pkgs_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

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

#helper functions
get_columns_by_prefix <- function(prefix) {
  grep(paste0('^', prefix), colnames(data), value = TRUE)
}

summarise_missing <- function(cols, label) {
  print(paste('The number of', label, 'are:', length(cols)))
  n_na <- sum(sapply(data[cols], function(x) any(is.na(x))))
  print(paste(label, 'number of columns with NA:', n_na))
}

convert_and_check_numeric <- function(cols) {
  for (col in cols) {
    if (!is.numeric(data[[col]])) {
      data[[col]] <<- as.numeric(data[[col]])
      if (any(is.na(data[[col]]))) {
        warning(paste('Some values in column', col,
                      'could not be converted to numeric and have been set to NA.'))
      }
    }
  }
  if (!all(sapply(data[cols], is.numeric))) {
    stop('Not all specified columns are numeric. Please check the data.')
  }
}

#features for metabolites and taxa
m_columns <- get_columns_by_prefix('m__')
t_columns <- get_columns_by_prefix('t__')

# check for missing values
summarise_missing(m_columns, 'metabolites')
summarise_missing(t_columns, 'taxa')

# check for response and stratify variables. 
print(paste('TRUE/FALSE: The response variable,', response_variable,
            'has missing values.', any(is.na(data[[response_variable]]))))

if (identical(stratify_variable, '')) {
  stratify_variable <- NULL
  print('Notice: The stratify variable is NULL.')
} else if (!is.null(stratify_variable)) {
  print(paste('TRUE/FALSE: The stratify variable,', stratify_variable,
              'has missing values.', any(is.na(data[[stratify_variable]]))))
}
#check for missing values in the entire dataset. NAs can throw errors in modeling. 
all_columns_to_check <- c(response_variable, m_columns, t_columns)
if (!is.null(stratify_variable)) {
  all_columns_to_check <- c(all_columns_to_check, stratify_variable)
}

if (any(sapply(data[all_columns_to_check], function(x) any(is.na(x))))) {
  stop('There are NA values in the dataset. Please check the data again before proceeding.')
} else {
  print('No missing values detected in the dataset. Analysis will proceed.')
}

# check whether classes are numeric. 
convert_and_check_numeric(c(m_columns, t_columns))
print('All metabolite and taxa columns are numeric.')

# Metrics and response variable type
if (type_of_analysis == 'binary') {
  source('model_functions/binary_functions.R')
  metrics <- list(
    accuracy = accuracy_calculation,
    kappa    = kappa_calculation,
    auroc   = auroc_calculation
  )
  print('For binary table this is the counts recorded:')
  print(table(data[[response_variable]]))
  print('Binary functions Accuracy, AUROC and Kappa have been loaded.')
  
  if (!is.factor(data[[response_variable]])) {
    data[[response_variable]] <- as.factor(data[[response_variable]])
    print(paste('The response variable', response_variable, 'has been converted to a factor.'))
  }
  
} else if (type_of_analysis == 'continuous') {
  source('model_functions/continuous_functions.R')
  metrics <- list(
    rmse = rmse_calculation,
    r2   = r2_calculation,
    mae  = mae_calculation
  )
  print('Continuous functions R^2, RMSE and MAE have been loaded.')
  
  if (!is.numeric(data[[response_variable]])) {
    data[[response_variable]] <- as.numeric(data[[response_variable]])
    print(paste('The response variable', response_variable, 'has been converted to numeric.'))
  }
  
} else {
  stop("Invalid type_of_analysis. Please enter 'binary' or 'continuous'.")
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

# Subset the testing data
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

