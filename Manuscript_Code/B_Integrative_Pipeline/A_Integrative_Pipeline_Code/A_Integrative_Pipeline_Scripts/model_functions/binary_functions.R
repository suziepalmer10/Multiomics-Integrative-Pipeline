#notes on inputs for these functions are working correctly
#Accuracy assumes that y_true and predictions are categorical with the same classes.
#AUROC assumes that y_true is a binary vector (e.g., 0/1) and predictions are probability
#scores or predicted probabilities.
#Kappa assumes that y_true and predictions are factors or categorical variables.

# CALCULATING ACCURACY VALUES
accuracy_calculation <- function(y_true, predictions) {
  predictions_ <- ifelse(predictions > 0.5, 1, 0)
  # Convert to numeric
  y_true <- as.numeric(as.character(y_true))
  predictions_ <- as.numeric(as.character(predictions_))
  accuracy <- mean(predictions_ == y_true)
  return(accuracy)
}

# CALCULATING AUC-ROC VALUES
auroc_calculation <- function(y_true, predictions) {
    roc_obj <- roc(y_true, predictions)
    auroc <- auc(roc_obj)
    return(auroc)
}

# CALCULATING KAPPA VALUES
kappa_calculation <- function(y_true, predictions) {
    # Ensure y_true and predictions are factors for kappa calculation
    predictions_ <- ifelse(predictions > 0.5, 1, 0)
    y_true <- as.numeric(y_true)
    predictions_ <- as.numeric(predictions_)
    pred_df <- data.frame(
      y_true_val = y_true,
      predictions_val = predictions_
    )
    kappa <- kappa2(pred_df, weight ='unweighted')
    return(kappa$value)
}


#load in lists for validation set
# Initialize lists for kappa
kappa_fold_metabolomics <- list()
kappa_fold_mss<- list()
kappa_fold_concat <- list()
# Initialize lists for accuracy
accuracy_fold_metabolomics <- list()
accuracy_fold_mss<- list()
accuracy_fold_concat <- list()
# Initialize lists for auroc
auroc_fold_metabolomics <- list()
auroc_fold_mss<- list()
auroc_fold_concat <- list()

#load in lists for testing set
#Initialize lists for kappa
test_kappa_fold_metabolomics <- list()
test_kappa_fold_mss<- list()
test_kappa_fold_concat <- list()
# Initialize lists for accuracy
test_accuracy_fold_metabolomics <- list()
test_accuracy_fold_mss<- list()
test_accuracy_fold_concat <- list()
#initialize lists for auroc
test_auroc_fold_metabolomics <- list()
test_auroc_fold_mss<- list()
test_auroc_fold_concat <- list()

#initialize lists for averaged stacked, weighted nnls and lasso/enet
#auroc
auroc_averaged_stacked <- list()
test_auroc_averaged_stacked <- list()
auroc_weighted_nnls <- list()
test_auroc_weighted_nnls <- list()
auroc_sparse_nnls <- list()
test_auroc_sparse_nnls <- list()
auroc_pls <- list()
test_auroc_pls <- list()

#accuracy
accuracy_averaged_stacked <- list()
test_accuracy_averaged_stacked <- list()
accuracy_weighted_nnls <- list()
test_accuracy_weighted_nnls <- list()
accuracy_sparse_nnls <- list()
test_accuracy_sparse_nnls <- list()
accuracy_pls <- list()
test_accuracy_pls <- list()

#kappa
kappa_averaged_stacked <- list()
test_kappa_averaged_stacked <- list()
kappa_weighted_nnls <- list()
test_kappa_weighted_nnls <- list()
kappa_sparse_nnls <- list()
test_kappa_sparse_nnls <- list()
kappa_pls <- list()
test_kappa_pls <- list()


