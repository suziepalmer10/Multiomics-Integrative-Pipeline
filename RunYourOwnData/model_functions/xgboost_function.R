model <- function(trainData, validationData, testingData, response_variable, response_type) {
  set.seed(123)
  
  # Prepare the data
  if (response_type == "binary") {
    trainData[[response_variable]] <- as.factor(trainData[[response_variable]])
    validationData[[response_variable]] <- as.factor(validationData[[response_variable]])
    testingData[[response_variable]] <- as.factor(testingData[[response_variable]])
    
    dtrain <- xgb.DMatrix(data = as.matrix(trainData %>% select(-all_of(response_variable))), 
                          label = as.numeric(trainData[[response_variable]]) - 1)
    dtest <- xgb.DMatrix(data = as.matrix(validationData %>% select(-all_of(response_variable))), 
                         label = as.numeric(validationData[[response_variable]]) - 1)
    dactual_test <- xgb.DMatrix(data = as.matrix(testingData %>% select(-all_of(response_variable))), 
                                label = as.numeric(testingData[[response_variable]]) - 1)
  } else if (response_type == "continuous") {
    trainData[[response_variable]] <- as.numeric(trainData[[response_variable]])
    validationData[[response_variable]] <- as.numeric(validationData[[response_variable]])
    testingData[[response_variable]] <- as.numeric(testingData[[response_variable]])
    
    dtrain <- xgb.DMatrix(data = as.matrix(trainData %>% select(-all_of(response_variable))), 
                          label = trainData[[response_variable]])
    dtest <- xgb.DMatrix(data = as.matrix(validationData %>% select(-all_of(response_variable))), 
                         label = validationData[[response_variable]])
    dactual_test <- xgb.DMatrix(data = as.matrix(testingData %>% select(-all_of(response_variable))), 
                                label = testingData[[response_variable]])
  } else {
    stop("Unsupported response type. Please use 'binary' or 'continuous'.")
  }
  
  # Define the metric
  metric <- if (response_type == "binary") "Accuracy" else "RMSE"
  
  # Set up the train control
  train_control <- trainControl(
    method = "cv", 
    number = 5,
    verboseIter = TRUE
  )

  # Stage 1: Tune core hyperparameters
  tune_grid_stage1 <- expand.grid(
    nrounds = c(100, 150, 200),
    max_depth = c(4, 6, 8, 10),
    eta = c(0.01, 0.1, 0.2),
    colsample_bytree = c(0.6, 0.8, 1),
    subsample = c(0.6, 0.8, 1),
    gamma = 0,
    min_child_weight = 1
  )
  
  model_stage1 <- caret::train(
    x = as.matrix(trainData %>% select(-all_of(response_variable))),
    y = if (response_type == "binary") as.factor(trainData[[response_variable]]) else trainData[[response_variable]],
    method = "xgbTree",
    metric = metric,
    tuneGrid = tune_grid_stage1,
    trControl = train_control
  )
  
  # Best parameters from Stage 1
  best_params_stage1 <- model_stage1$bestTune
  
  # Stage 2: Fine-tune additional parameters
  tune_grid_stage2 <- expand.grid(
    nrounds = best_params_stage1$nrounds,
    max_depth = best_params_stage1$max_depth,
    eta = best_params_stage1$eta,
    colsample_bytree = best_params_stage1$colsample_bytree,
    subsample = best_params_stage1$subsample,
    gamma = c(0, 1, 5),
    min_child_weight = c(1, 5, 10)
  )
  
  model_stage2 <- caret::train(
    x = as.matrix(trainData %>% select(-all_of(response_variable))),
    y = if (response_type == "binary") as.factor(trainData[[response_variable]]) else trainData[[response_variable]],
    method = "xgbTree",
    metric = metric,
    tuneGrid = tune_grid_stage2,
    trControl = train_control
  )
  
  # Make predictions on validationData
  predictions <- if (response_type == "binary") {
    predict(model_stage2, newdata = as.matrix(validationData %>% select(-all_of(response_variable))), type = "prob")[, 2]
  } else {
    predict(model_stage2, newdata = as.matrix(validationData %>% select(-all_of(response_variable))))
  }
  
  predictions_response <- if (response_type == "binary") {
    as.numeric(validationData[[response_variable]] == levels(validationData[[response_variable]])[2])
  } else {
    validationData[[response_variable]]
  }
  
  # Make predictions on testingData
  test_predictions <- if (response_type == "binary") {
    predict(model_stage2, newdata = as.matrix(testingData %>% select(-all_of(response_variable))), type = "prob")[, 2]
  } else {
    predict(model_stage2, newdata = as.matrix(testingData %>% select(-all_of(response_variable))))
  }
  
  test_predictions_response <- if (response_type == "binary") {
    as.numeric(testingData[[response_variable]] == levels(testingData[[response_variable]])[2])
  } else {
    testingData[[response_variable]]
  }
  
  return(list(predictions = predictions, 
              predictions_response = predictions_response, 
              test_predictions = test_predictions, 
              test_predictions_response = test_predictions_response, 
              model = model_stage2))
}

