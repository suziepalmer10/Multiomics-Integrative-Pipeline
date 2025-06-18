#This file is specifically for evaluating performance metrics for continuous data

# Combine the lists into a list of lists
kappa_train_results <- list(kappa_fold_metabolomics, kappa_fold_mss, kappa_fold_concat, kappa_averaged_stacked, kappa_weighted_nnls, kappa_sparse_nnls, kappa_pls)
auroc_train_results <- list(auroc_fold_metabolomics, auroc_fold_mss, auroc_fold_concat, auroc_averaged_stacked, auroc_weighted_nnls, auroc_sparse_nnls, auroc_pls)
accuracy_train_results <- list(accuracy_fold_metabolomics, accuracy_fold_mss, accuracy_fold_concat, accuracy_averaged_stacked, accuracy_weighted_nnls, accuracy_sparse_nnls, accuracy_pls)
kappa_test_results <- list(test_kappa_fold_metabolomics, test_kappa_fold_mss, test_kappa_fold_concat, test_kappa_averaged_stacked, test_kappa_weighted_nnls, test_kappa_sparse_nnls, test_kappa_pls)
auroc_test_results <- list(test_auroc_fold_metabolomics, test_auroc_fold_mss, test_auroc_fold_concat, test_auroc_averaged_stacked, test_auroc_weighted_nnls, test_auroc_sparse_nnls, test_auroc_pls)
accuracy_test_results <- list(test_accuracy_fold_metabolomics, test_accuracy_fold_mss, test_accuracy_fold_concat, test_accuracy_averaged_stacked, test_accuracy_weighted_nnls, test_accuracy_sparse_nnls, test_accuracy_pls)

# Define the row labels
row_labels <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS")

# Get the dataframes for each metric
kappa_train_df <- stats_results(kappa_train_results, row_labels, "Kappa Train")
auroc_train_df <- stats_results(auroc_train_results, row_labels, "auroc Train")
accuracy_train_df <- stats_results(accuracy_train_results, row_labels, "Accuracy Train")
kappa_test_df <- stats_results(kappa_test_results, row_labels, "Kappa Test")
auroc_test_df <- stats_results(auroc_test_results, row_labels, "auroc Test")
accuracy_test_df <- stats_results(accuracy_test_results, row_labels, "Accuracy Test")

# Merge all dataframes by 'Category'
all_results_df <- reduce(list(kappa_train_df, auroc_train_df, accuracy_train_df, kappa_test_df, auroc_test_df, accuracy_test_df), function(x, y) {
  full_join(x, y, by = "Category")
})

# Define the order for the facets and x-axis categories
facet_order <- c("Kappa Train ", "Accuracy Train ", "AUROC Train ", "Kappa Test ", "Accuracy Test ", "AUROC Test ") # Replace with your actual metric names in the desired order
x_order <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked", "Weighted NNLS", "Lasso Stacked", "PLS") # Replace with your actual category names in the desired order
