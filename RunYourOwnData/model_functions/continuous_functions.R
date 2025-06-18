# CALCULATING R^2 VALUES
r2_calculation <- function(y_true, predictions) {
  TSS <- sum((y_true - mean(y_true))^2)
  RSS <- sum((y_true - predictions)^2)
  R2 <- 1 - (RSS / TSS)
  return(R2)
}

# CALCULATING RMSE VALUES
rmse_calculation <- function(y_true, predictions) {
  MSE <- mean((y_true - predictions)^2)
  RMSE <- sqrt(MSE)
  return(RMSE)
}

# CALCULATING MAE VALUES
mae_calculation <- function(y_true, predictions) {
  MAE <- mean(abs(y_true - predictions))
  return(MAE)
}

#load in lists for validation set
# Initialize lists for rmse
rmse_fold_metabolomics <- list()
rmse_fold_mss<- list()
rmse_fold_concat <- list()
# Initialize lists for r2
r2_fold_metabolomics <- list()
r2_fold_mss<- list()
r2_fold_concat <- list()
# Initialize lists for mae
mae_fold_metabolomics <- list()
mae_fold_mss<- list()
mae_fold_concat <- list()

#load in lists for testing set
#Initialize lists for rmse
test_rmse_fold_metabolomics <- list()
test_rmse_fold_mss<- list()
test_rmse_fold_concat <- list()
# Initialize lists for r2
test_r2_fold_metabolomics <- list()
test_r2_fold_mss<- list()
test_r2_fold_concat <- list()
#initialize lists for mae
test_mae_fold_metabolomics <- list()
test_mae_fold_mss<- list()
test_mae_fold_concat <- list()

#initialize lists for averaged stacked, weighted nnls and lasso/enet
#rmse
rmse_averaged_stacked <- list()
test_rmse_averaged_stacked <- list()
rmse_weighted_nnls <- list()
test_rmse_weighted_nnls <- list()
rmse_sparse_nnls <- list()
test_rmse_sparse_nnls <- list()
rmse_pls <- list()
test_rmse_pls <-list()
#r2
r2_averaged_stacked <- list()
test_r2_averaged_stacked <- list()
r2_weighted_nnls <- list()
test_r2_weighted_nnls <- list()
r2_sparse_nnls <- list()
test_r2_sparse_nnls <- list()
r2_pls <- list()
test_r2_pls <- list()
#mae
mae_averaged_stacked <- list()
test_mae_averaged_stacked <- list()
mae_weighted_nnls <- list()
test_mae_weighted_nnls <- list()
mae_sparse_nnls <- list()
test_mae_sparse_nnls <- list()
mae_pls <- list()
test_mae_pls <- list()













