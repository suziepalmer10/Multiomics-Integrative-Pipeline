#dataset: Erawijantari
#purpose: data preparation
library("data.table")
library(dplyr)
library(tidyverse)
library(tidyr)
library(limma)
library(openxlsx)
library(writexl)

setwd("/home2/s180020/Desktop/Franzosa_integrated/OriginalDataProcessing")

#SECTION 1: COMBINE, SCALE AND FILTER DATA. 
metadata_cols <- c("Sample", "Study.Group", 'disease_status', 'Fecal.Calprotectin', 'Age')
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
counts <- table(metadata$Study.Group)
print(counts)

#metadata to keep
#add an additional column, which converts our data into binary Control or Disease. 
metadata <- metadata %>%
  mutate(disease_status = ifelse(Study.Group %in% c("UC", "CD"), "Disease","Control"))

#add m__ before each metabolite. 
#this will allow us to parse and automate the data downstream
metabolites <- metabolites %>%
  rename_with(~ paste0("m__", .), -Sample)
# Apply log2 transformation only to numeric columns, excluding "Sample"
metabolites <- metabolites %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))

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
metadata$Fecal.Calprotectin <- log2(metadata$Fecal.Calprotectin)
summary(metadata$Fecal.Calprotectin)


metadata_to_merge <- metadata[, ..metadata_cols]
#the code below counts the number of NAs for the Fecal.Calprotectin groups, 
#since there is 67 missing values for these. 
metadata[is.na(Fecal.Calprotectin), .N, by = Study.Group]
#NA for age? 
metadata[is.na(Age)]
#No missing values for Age

#merge metadata and metabolites
metabolites_with_labels <- merge(metadata_to_merge, metabolites, by = "Sample")
franzosa_processed_data <- merge(metabolites_with_labels, filtered_mss, by = "Sample")

##THIS SECTION IS TO GENERATE non-processed datasets (no differential expression - as controls).
franzosa_convsdis_age_node <- franzosa_processed_data
franzosa_convsdis_fp_node <- franzosa_convsdis_age_node  %>% filter(!is.na(Fecal.Calprotectin))
franzosa_convsUC_age_node <- franzosa_convsdis_age_node %>% filter(Study.Group != "CD")
franzosa_convsUC_fp_node <- franzosa_convsdis_fp_node %>% filter(Study.Group != "CD")
franzosa_convsCD_age_node <- franzosa_convsdis_age_node %>% filter(Study.Group != "UC")
franzosa_convsCD_fp_node <- franzosa_convsdis_fp_node %>% filter(Study.Group != "UC")

write.csv(franzosa_convsdis_age_node, file = "nofilter_franzosa_ControlvsDisease_Age_11122024.csv", row.names = FALSE) 
write.csv(franzosa_convsdis_fp_node, file = "nofilter_franzosa_ControlvsDisease_Fecal_Calprotectin_11122024.csv", row.names = FALSE) 
write.csv(franzosa_convsUC_age_node, file = "nofilter_franzosa_ControlvsUC_Age_11122024.csv", row.names = FALSE) 
write.csv(franzosa_convsUC_fp_node, file = "nofilter_franzosa_ControlvsUC_Fecal_Calprotectin_11122024.csv", row.names = FALSE) 
write.csv(franzosa_convsCD_age_node, file = "nofilter_franzosa_ControlvsCD_Age_11122024.csv", row.names = FALSE) 
write.csv(franzosa_convsCD_fp_node, file = "nofilter_franzosa_ControlvsCD_Fecal_Calprotectin_11122024.csv", row.names = FALSE) 



#SECTION 2: DIFFERENTIAL EXPRESSION ANALYSIS FOR FEATURE REDUCTION. 
#will perform separately for each dataset. 

#PERFORM LIMMA FOR METABOLITES
m_columns <- grep("^m__", colnames(franzosa_processed_data), value = TRUE)
metabolite_data <- franzosa_processed_data %>% select(m_columns)
rownames(metabolite_data) <- franzosa_processed_data$Sample
metabolite_data <- as.matrix(metabolite_data)
metabolite_data_t <- t(metabolite_data)
metadata_metab <- franzosa_processed_data %>% 
  select(Study.Group) %>% 
  as.data.frame()
# Set row names using the Sample column from franzosa_processed_data
rownames(metadata_metab) <- franzosa_processed_data$Sample

#use metadata, design and contrasts from above for this problem. 
#contrasts fit
# Fit the linear model
design <- model.matrix(~ Study.Group-1, data = metadata_metab)
# Define contrasts
contrasts <- makeContrasts(
  Control_vs_UC = Study.GroupControl - Study.GroupUC,
  Control_vs_CD = Study.GroupControl - Study.GroupCD,
  levels = design
)
fit_m <- lmFit(metabolite_data_t, design)
# Apply contrasts
fit_contrasts_m <- contrasts.fit(fit_m, contrasts)
# Apply empirical Bayes moderation
fit_contrasts_m <- eBayes(fit_contrasts_m)

# # Define a function to extract and filter results based on adjusted p-value and fold change
extract_significant_taxa <- function(fit_contrasts, contrast_name, p_value_threshold = 0.05, fold_change_threshold = 2.0) {
  # Extract results for the given contrast
  results <- topTable(fit_contrasts, coef = contrast_name, adjust = "BH", p.value = p_value_threshold, number = Inf)
  # Filter results based on adjusted p-value and fold change
  significant_hits <- results[
    results$adj.P.Val < p_value_threshold &
      (results$logFC > fold_change_threshold | results$logFC < -fold_change_threshold),
  ]
  return(significant_hits)
}

#List of contrasts to process
contrast_names <- c(
  "Control_vs_UC",
  "Control_vs_CD"
)

# Apply the function to each contrast and store the results in a list
results_list_m <- lapply(contrast_names, function(contrast) {
  extract_significant_taxa(fit_contrasts_m, contrast)
})
# Optionally, assign names to the list elements for easier access
names(results_list_m) <- contrast_names

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
saveWorkbook(wb, file = "Franzosa_Differential_Abundance_Results_with_Metabolites_2.0fc_0.05pval_10112024.xlsx", overwrite = TRUE)

#these will be uncoupled in three ways based on ControlvDisease, ControlvUC and ControlvCD. 
#reduce features to significant only. We are going to use these to generate some models.
significant_metabolites <- lapply(results_list_m, function(df) {
  rownames_list <- rownames(df)
  return(rownames(df))
})
#this will be used for the metabolomics data. 
diff_expressed_metabolites_all <- unique(unlist(significant_metabolites, recursive = FALSE))
diff_expressed_controlvUC <- unlist(rownames(results_list_m[[1]]))
diff_expressed_controlvCD <- unlist(rownames(results_list_m[[2]]))


#PERFORM WILCOXON FOR TAXA
t_columns <- grep("^t__", colnames(franzosa_processed_data), value = TRUE)
w_mss <- franzosa_processed_data %>% select('Study.Group', t_columns)
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
treatment_groups <- c("UC", "CD")

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
diff_expressed_taxa_ControlvsUC<- results %>%
  filter(Comparison == "Control vs UC")
diff_expressed_taxa_ControlvsUC <- unlist(diff_expressed_taxa_ControlvsUC$Taxon)
diff_expressed_taxa_ControlvsCD<- results %>%
  filter(Comparison == "Control vs CD")
diff_expressed_taxa_ControlvsCD <- unlist(diff_expressed_taxa_ControlvsCD$Taxon)

split_data <- split(results, results$Comparison)
names(split_data) <- gsub("/", "_", names(split_data))
# Specify the output Excel file path
output_file <- "franzosa_Differential_Abundance_Results_with_Taxa_0.05pval_11122024.xlsx"
# Write each dataframe to a separate sheet in the Excel file
write_xlsx(split_data, path = output_file)

#SECTION 3: COMBINE ABX METAB AND TAXA - COMBINE DIET METAB AND TAXA

#subset for disease status (i.e., control vs disease)
diff_expressed_data_list_all <- c(metadata_cols, diff_expressed_metabolites_all, diff_expressed_taxa_all)
diff_expressed_data_all <- franzosa_processed_data[, ..diff_expressed_data_list_all]
diff_expressed_data_all_fp <- diff_expressed_data_all %>%
  filter(!is.na(Fecal.Calprotectin))
write.csv(diff_expressed_data_all, file = "franzosa_ControlvsDisease_Age_11122024.csv", row.names = FALSE) 
write.csv(diff_expressed_data_all_fp, file = "franzosa_ControlvsDisease_Fecal_Calprotectin_11122024.csv", row.names = FALSE) 

#subset for control vs UC
diff_expressed_data_list_UC <- c(metadata_cols, diff_expressed_controlvUC, diff_expressed_taxa_ControlvsUC)
diff_expressed_data_UC <- franzosa_processed_data[, ..diff_expressed_data_list_UC]
diff_expressed_data_UC<- diff_expressed_data_UC %>%
  filter(Study.Group != "CD")
diff_expressed_data_UC_fp<- diff_expressed_data_UC %>%
  filter(!is.na(Fecal.Calprotectin))
write.csv(diff_expressed_data_UC, file = "franzosa_ControlvsUC_Age_11122024.csv", row.names = FALSE) 
write.csv(diff_expressed_data_UC_fp, file = "franzosa_ControlvsUC_Fecal_Calprotectin_11122024.csv", row.names = FALSE) 

#subset for control vs CD
diff_expressed_data_list_CD <- c(metadata_cols, diff_expressed_controlvCD, diff_expressed_taxa_ControlvsCD)
diff_expressed_data_CD <- franzosa_processed_data[, ..diff_expressed_data_list_CD]
diff_expressed_data_CD<- diff_expressed_data_CD %>%
  filter(Study.Group != "UC")
diff_expressed_data_CD_fp<- diff_expressed_data_CD %>%
  filter(!is.na(Fecal.Calprotectin))
write.csv(diff_expressed_data_CD, file = "franzosa_ControlvsCD_Age_11122024.csv", row.names = FALSE) 
write.csv(diff_expressed_data_CD_fp, file = "franzosa_ControlvsCD_Fecal_Calprotectin_11122024.csv", row.names = FALSE) 
