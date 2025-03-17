model <- function(trainData, testData, actual_testing_Data, response_variable, response_type){
  set.seed(123)
  
  # Define the tuning grid for mtry
  tune_grid <- expand.grid(mtry = seq(2, sqrt(ncol(trainData) - 1), by = 10))
  
  # Determine the method and family based on response type
  if (response_type == "binary") {
    method <- "rf"
    metric <- "Accuracy"
  } else if (response_type == "continuous") {
    method <- "rf"
    metric <- "RMSE"
  } else {
    stop("Unsupported response type. Please use 'binary' or 'continuous'.")
  }
  
  # Train the model
  model <- caret::train(
    x = trainData[, -which(names(trainData) == response_variable)],
    y = trainData[[response_variable]],
    method = method,
    metric = metric,
    tuneGrid = tune_grid,
    trControl = trainControl(method = "cv", number = 5)
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
