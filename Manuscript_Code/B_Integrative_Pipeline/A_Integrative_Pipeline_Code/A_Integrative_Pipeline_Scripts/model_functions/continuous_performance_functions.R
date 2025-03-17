#This file is specifically for evaluating performance metrics for continuous data

# Combine the lists into a list of lists
rmse_train_results <- list(rmse_fold_metabolomics, rmse_fold_mss, rmse_fold_concat, rmse_averaged_stacked, rmse_weighted_nnls, rmse_sparse_nnls, rmse_pls)
r2_train_results <- list(r2_fold_metabolomics, r2_fold_mss, r2_fold_concat, r2_averaged_stacked, r2_weighted_nnls, r2_sparse_nnls, r2_pls)
mae_train_results <- list(mae_fold_metabolomics, mae_fold_mss, mae_fold_concat, mae_averaged_stacked, mae_weighted_nnls, mae_sparse_nnls, mae_pls)
rmse_test_results <- list(test_rmse_fold_metabolomics, test_rmse_fold_mss, test_rmse_fold_concat, test_rmse_averaged_stacked, test_rmse_weighted_nnls, test_rmse_sparse_nnls, test_rmse_pls)
r2_test_results <- list(test_r2_fold_metabolomics, test_r2_fold_mss, test_r2_fold_concat, test_r2_averaged_stacked, test_r2_weighted_nnls, test_r2_sparse_nnls, test_r2_pls)
mae_test_results <- list(test_mae_fold_metabolomics, test_mae_fold_mss, test_mae_fold_concat, test_mae_averaged_stacked, test_mae_weighted_nnls, test_mae_sparse_nnls, test_mae_pls)

# Define the row labels
row_labels <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS")

# Get the dataframes for each metric
rmse_train_df <- stats_results(rmse_train_results, row_labels, "RMSE Train")
r2_train_df <- stats_results(r2_train_results, row_labels, "R2 Train")
mae_train_df <- stats_results(mae_train_results, row_labels, "MAE Train")
rmse_test_df <- stats_results(rmse_test_results, row_labels, "RMSE Test")
r2_test_df <- stats_results(r2_test_results, row_labels, "R2 Test")
mae_test_df <- stats_results(mae_test_results, row_labels, "MAE Test")

# Merge all dataframes by 'Category'
all_results_df <- reduce(list(rmse_train_df, r2_train_df, mae_train_df, rmse_test_df, r2_test_df, mae_test_df), function(x, y) {
  full_join(x, y, by = "Category")
})

# Define the order for the facets and x-axis categories
facet_order <- c("RMSE Train ", "MAE Train ", "R2 Train ", "RMSE Test ", "MAE Test ", "R2 Test ") # Replace with your actual metric names in the desired order
x_order <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS") # Replace with your actual category names in the desired order
