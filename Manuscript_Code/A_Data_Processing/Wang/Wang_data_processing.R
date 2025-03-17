#dataset: Wang
#purpose: data preparation
library("data.table")
library(dplyr)
library(tidyverse)
library(tidyr)
library(limma)
library(openxlsx)
library(writexl)

setwd("/home2/s180020/Desktop/Wang_integrated/OriginalDataProcessing")

#SECTION 1: COMBINE, SCALE AND FILTER DATA. 
metadata_cols <- c("Sample", "Study.Group", "Age", "Gender", "BMI", "Creatinine", "Urea", "eGFR")
##load in datasets
#load in genera relative abundance
genera <- fread("genera.tsv", header = TRUE)
#load in species relative abundance
species <- fread("species.tsv", header=TRUE)
#load in metabolites
metabolites <- fread("mtb.tsv", header=TRUE) 
#metabolite annotations
mapped_metabolites <- fread("mtb.map.tsv", header=TRUE) 
#metadata
metadata <- fread("metadata.tsv", header=TRUE)



#check number of each condition
table(metadata$Study.Group)
table(metadata$Gender)
summary(metadata$Age)
summary(metadata$BMI) #7 NA's
summary(metadata$Creatinine) #6 NA's
summary(metadata$Urea) #6 NA's
summary(metadata$eGFR) #6 NA's

#metadata to keep
#add m__ before each metabolite. 
#this will allow us to parse and automate the data downstream
metabolites <- metabolites %>%
  rename_with(~ paste0("m__", .), -Sample)
# # Apply log2 transformation only to numeric columns, excluding "Sample"
# metabolites <- metabolites %>%
#   mutate(across(where(is.numeric), ~ log2(. + 1)))

#merge species and genus by Sample
mss <- merge(genera, species, by = "Sample")
# Calculate sums of numeric columns
sums <- mss %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
# Select columns with sums >= 0.001
filtered_mss <- mss %>% select(where(~ !is.numeric(.) || sum(., na.rm = TRUE) >= 0.001))
filtered_sums <- filtered_mss %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))

#add t__ before all taxa to make easier to parse
filtered_mss<- filtered_mss %>%
  rename_with(~ paste0("t__", .), -Sample)

#log transform Fecal.calprotectin
metadata$Creatinine <- log2(metadata$Creatinine)
summary(metadata$Creatinine)

metadata_to_merge <- metadata[, ..metadata_cols]
#the code below counts the number of NAs for the Fecal.Calprotectin groups, 
#checked above for NAS

#merge metadata and metabolites
metabolites_with_labels <- merge(metadata_to_merge, metabolites, by = "Sample")
wang_processed_data <- merge(metabolites_with_labels, filtered_mss, by = "Sample")

##THIS SECTION IS TO GENERATE non-processed datasets (no differential expression - as controls).
wang_studygroup_age_gender <- wang_processed_data
wang_bmi <- wang_studygroup_age_gender  %>% filter(!is.na(BMI))
wang_egfr_creatinine_urea <- wang_studygroup_age_gender %>%
  filter(!is.na(eGFR), !is.na(Creatinine), !is.na(Urea))

write.csv(wang_studygroup_age_gender, file = "nofilter_wang_studygroup_age_gender_12092024.csv", row.names = FALSE) 
write.csv(wang_bmi, file = "nofilter_wang_bmi_12092024.csv", row.names = FALSE)
write.csv(wang_egfr_creatinine_urea, file = "nofilter_wang_egfr_creatinine_urea.csv", row.names = FALSE)

#SECTION 2: DIFFERENTIAL EXPRESSION ANALYSIS FOR FEATURE REDUCTION. 
#will perform separately for each dataset. 

#PERFORM LIMMA FOR METABOLITES
m_columns <- grep("^m__", colnames(wang_processed_data), value = TRUE)
metabolite_data <- wang_processed_data %>% select(m_columns)
rownames(metabolite_data) <- wang_processed_data$Sample
metabolite_data <- as.matrix(metabolite_data)
metabolite_data_t <- t(metabolite_data)
metadata_metab <- wang_processed_data %>% 
  select(Study.Group) %>% 
  as.data.frame()
# Set row names using the Sample column from franzosa_processed_data
rownames(metadata_metab) <- wang_processed_data$Sample

#use metadata, design and contrasts from above for this problem. 
#contrasts fit
# Fit the linear model
design <- model.matrix(~ Study.Group-1, data = metadata_metab)
# Define contrasts
contrasts <- makeContrasts(
  Control_vs_ESRD = Study.GroupControl - Study.GroupESRD,
  levels = design
)
fit_m <- lmFit(metabolite_data_t, design)
# Apply contrasts
fit_contrasts_m <- contrasts.fit(fit_m, contrasts)
# Apply empirical Bayes moderation
fit_contrasts_m <- eBayes(fit_contrasts_m)

# Define a function to extract and filter results based on adjusted p-value and fold change
extract_significant_taxa <- function(fit_contrasts, contrast_name, p_value_threshold = 0.05, fold_change_threshold = 2.0, filter_results = TRUE) {
  # Extract results for the given contrast
  results <- topTable(fit_contrasts, coef = contrast_name, adjust = "BH", p.value = 1, number = Inf)
  
  if (filter_results) {
    # Filter results based on adjusted p-value and fold change
    significant_hits <- results[
      #results$P.Value <p_value_threshold,
      results$adj.P.Val < p_value_threshold, #&
        #(results$logFC > fold_change_threshold | results$logFC < -fold_change_threshold),
    ]
    return(significant_hits)
  } else {
    # Return all results without filtering
    return(results)
  }
}


#List of contrasts to process
contrast_names <- c(
  "Control_vs_ESRD"
)

# Apply the function to each contrast and store the results in a list
results_list_m <- lapply(contrast_names, function(contrast) {
  extract_significant_taxa(fit_contrasts_m, contrast)
})


# Optionally, assign names to the list elements for easier access
names(results_list_m) <- contrast_names

df <- results_list_m[[1]]

library(ggplot2)

# Example function to create a volcano plot
create_volcano_plot <- function(results) {
  # Ensure the necessary columns exist
  if (!all(c("logFC", "adj.P.Val") %in% colnames(results))) {
    stop("The results dataframe must contain 'logFC' and 'adj.P.Val' columns.")
  }
  
  # Add a column for point color based on the p-value
  results$Color <- ifelse(-log10(results$adj.P.Val) > 2, "blue", "black")
  
  # Create the volcano plot
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Color)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_identity() +  # Use the colors directly from the 'Color' column
    theme_minimal() +
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value"
    ) +
    theme(
      legend.position = "none",  # Hide legend as we define colors explicitly
      plot.title = element_text(hjust = 0.5)
    ) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray")
  
  return(volcano_plot)
}

# Example usage (assuming 'results' is your dataframe from the function)
# results <- extract_significant_taxa(...)
volcano_plot <- create_volcano_plot(df)
print(volcano_plot)

df


# Create a new Excel workbook
wb <- createWorkbook()
# Loop through the list and add each data frame to a new sheet
for (i in seq_along(results_list_m)) {
  # Convert row names metabolites to a column
  df <- results_list_m[[i]]
  df <- cbind(Metabolites = rownames(df), df)
  rownames(df) <- NULL  # Remove the original row names
  # Add a worksheet with the name of the contrast
  addWorksheet(wb, sheetName = names(results_list_m)[i])
  # Write the data frame with taxa names to the corresponding sheet
  writeData(wb, sheet = names(results_list_m)[i], df)
}

# Save the workbook to a file
saveWorkbook(wb, file = "Wang_Differential_Abundance_Results_with_Metabolites_0.05adjpval_12102024.xlsx", overwrite = TRUE)

#reduce features to significant only. We are going to use these to generate some models.
significant_metabolites <- lapply(results_list_m, function(df) {
  rownames_list <- rownames(df)
  return(rownames(df))
})
#this will be used for the metabolomics data. 
diff_expressed_metabolites_all <- unique(unlist(significant_metabolites, recursive = FALSE))



#PERFORM WILCOXON FOR TAXA
t_columns <- grep("^t__", colnames(wang_processed_data), value = TRUE)
w_mss <- wang_processed_data %>% select('Study.Group', t_columns)
# Convert mss to a dataframe if it's a tibble
w_mss <- as.data.frame(w_mss)

# Initialize results data frame
results <- data.frame(Taxon = character(), Comparison = character(), p_value = numeric())
# Function to calculate fold change and Wilcoxon test
perform_comparison <- function(group1, group2, taxon) {
  # Use pull() to ensure a numeric vector is extracted
  group1_data <- as.numeric(w_mss %>% filter(Study.Group == group1) %>% pull(!!sym(taxon)))
  group2_data <- as.numeric(w_mss %>% filter(Study.Group == group2) %>% pull(!!sym(taxon)))
  # Calculate median values for each group
  median_group1 <- median(group1_data, na.rm = TRUE)
  median_group2 <- median(group2_data, na.rm = TRUE)
  # Perform Wilcoxon test
  test_result <- wilcox.test(group1_data, group2_data, exact=FALSE)
  # Return results
  return(data.frame(
    Taxon = taxon,
    Comparison = paste(group1, "vs", group2),
    p_value = test_result$p.value
  ))
}
# Define control and comparison groups
control_group <- "Control"
treatment_groups <- 'ESRD'

# Loop through each taxon and perform comparisons
for (taxon in colnames(w_mss)[2:ncol(w_mss)]) {  # Adjust starting column as needed
  # Control vs. treatment groups
  for (treatment in treatment_groups) {
    results <- rbind(results, perform_comparison(control_group, treatment, taxon))
  }
}

#adjust the p-value for multiple testing. 
results$adjusted_p <- p.adjust(results$p_value, method = "BH")

#filter out adjusted p-values that are >=0.05 or NA
results <- results[!is.na(results$adjusted_p) & results$adjusted_p <= 0.05, ]

#differentially expressed taxa list
diff_expressed_taxa_all <- unlist(unique(results$Taxon))


split_data <- split(results, results$Comparison)
names(split_data) <- gsub("/", "_", names(split_data))
# Specify the output Excel file path
output_file <- "wang_Differential_Abundance_Results_with_Taxa_0.05pval_12092024.xlsx"
# Write each dataframe to a separate sheet in the Excel file
write_xlsx(split_data, path = output_file)

#SECTION 3: COMBINE ABX METAB AND TAXA - COMBINE DIET METAB AND TAXA



#subset for disease status (i.e., control vs disease)
diff_expressed_data_list_all <- c(metadata_cols, diff_expressed_metabolites_all, diff_expressed_taxa_all)
diff_expressed_data_all <- wang_processed_data[, ..diff_expressed_data_list_all]

##THIS SECTION IS TO GENERATE non-processed datasets (no differential expression - as controls).
filtered_wang_studygroup_age_gender <- diff_expressed_data_all
filtered_wang_bmi <- filtered_wang_studygroup_age_gender  %>% filter(!is.na(BMI))
filtered_wang_egfr_creatinine_urea <- filtered_wang_studygroup_age_gender %>%
  filter(!is.na(eGFR), !is.na(Creatinine), !is.na(Urea))

write.csv(filtered_wang_studygroup_age_gender, file = "wang_studygroup_age_gender_12092024.csv", row.names = FALSE) 
write.csv(filtered_wang_bmi, file = "wang_bmi_12092024.csv", row.names = FALSE)
write.csv(filtered_wang_egfr_creatinine_urea, file = "wang_egfr_creatinine_urea.csv", row.names = FALSE)