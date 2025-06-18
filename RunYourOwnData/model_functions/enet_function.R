model <- function(trainData, validationData, testingData, response_variable, response_type) {
  set.seed(123)
  
  # Convert tibbles to data frames if needed
  trainData <- as.data.frame(trainData)
  validationData <- as.data.frame(validationData)
  testingData <- as.data.frame(testingData)
  
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
  
  # Make predictions on validationData
  if (response_type == "binary") {
    predictions <- predict(model, newdata = validationData[, -which(names(validationData) == response_variable)], type = "prob")[, 2]
    predictions_response <- as.numeric(validationData[[response_variable]] == levels(validationData[[response_variable]])[2])
  } else {
    predictions <- predict(model, newdata = validationData[, -which(names(validationData) == response_variable)])
    predictions_response <- validationData[[response_variable]]
  }
  
  # Make predictions on testingData
  if (response_type == "binary") {
    test_predictions <- predict(model, newdata = testingData[, -which(names(testingData) == response_variable)], type = "prob")[, 2]
    test_predictions_response <- as.numeric(testingData[[response_variable]] == levels(testingData[[response_variable]])[2])
  } else {
    test_predictions <- predict(model, newdata = testingData[, -which(names(testingData) == response_variable)])
    test_predictions_response <- testingData[[response_variable]]
  }
  
  return(list(predictions = predictions, 
              predictions_response = predictions_response, 
              test_predictions = test_predictions, 
              test_predictions_response = test_predictions_response, 
              model = model))
}


