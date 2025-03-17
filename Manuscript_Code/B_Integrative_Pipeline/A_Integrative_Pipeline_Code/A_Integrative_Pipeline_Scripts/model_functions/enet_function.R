model <- function(trainData, testData, actual_testing_Data, response_variable, response_type) {
  set.seed(123)
  
  # Convert tibbles to data frames if needed
  trainData <- as.data.frame(trainData)
  testData <- as.data.frame(testData)
  actual_testing_Data <- as.data.frame(actual_testing_Data)
  
 # Create a tuning grid
 #tune_grid <- expand.grid(
 #   alpha = seq(0, 1, by = 0.1),
 #   lambda = seq(0.01, 1, length = 10)
 # )
tune_grid <- expand.grid(
    alpha = seq(0, 1, by = 0.05), # More granularity in alpha
    lambda = seq(0.001, 10, length = 20) # Wider range and finer granularity for lambda
  )
  
  # Determine the method and family based on response type
  if (response_type == "binary") {
    method <- "glmnet"
    family <- "binomial"
  } else if (response_type == "continuous") {
    method <- "glmnet"
    family <- "gaussian"
  } else {
    stop("Unsupported response type. Please use 'binary' or 'continuous'.")
  }
  
  # Train the model
  model <- caret::train(
    x = trainData[, -which(names(trainData) == response_variable)],
    y = trainData[[response_variable]],
    method = method,
    family = family,
    tuneGrid = tune_grid
  )
  
  # Make predictions on testData
  if (response_type == "binary") {
    predictions <- predict(model, newdata = testData[, -which(names(testData) == response_variable)], type = "prob")[, 2]
    predictions_response <- as.numeric(testData[[response_variable]] == levels(testData[[response_variable]])[2])
  } else {
    predictions <- predict(model, newdata = testData[, -which(names(testData) == response_variable)])
    predictions_response <- testData[[response_variable]]
  }
  
  # Make predictions on actual_testing_Data
  if (response_type == "binary") {
    test_predictions <- predict(model, newdata = actual_testing_Data[, -which(names(actual_testing_Data) == response_variable)], type = "prob")[, 2]
    test_predictions_response <- as.numeric(actual_testing_Data[[response_variable]] == levels(actual_testing_Data[[response_variable]])[2])
  } else {
    test_predictions <- predict(model, newdata = actual_testing_Data[, -which(names(actual_testing_Data) == response_variable)])
    test_predictions_response <- actual_testing_Data[[response_variable]]
  }
  
  return(list(predictions = predictions, 
              predictions_response = predictions_response, 
              test_predictions = test_predictions, 
              test_predictions_response = test_predictions_response, 
              model = model))
}


