#dataset: Yachida
#purpose: data preparation
library("data.table")
library(dplyr)
library(tidyverse)
library(tidyr)
library(limma)
library(openxlsx)
library(writexl)

setwd("/home2/s180020/Desktop/Yachida_integrated/OriginalDataProcessing")

#SECTION 1: COMBINE, SCALE AND FILTER DATA. 
metadata_cols <- c("Sample", "Study.Group", 'disease_status', 'Age', 'Gender', 'Alcohol', 'BMI', 'Brinkman Index', 'healthy_vs_early')
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

#metadata for Control vs Cancer
#removes HS (history of colon surgery) and MP (multiple polyps)
metadata <- metadata %>%
  mutate(
    disease_status = ifelse(Study.Group %in% c("HS", "MP"), NA, 
                            ifelse(Study.Group %in% c("Stage_0", "Stage_I_II", "Stage_III_IV"), "Cancer", "Healthy"))
  )
#Control_vs_MP_Stage_0
metadata <- metadata %>%
  mutate(
    healthy_vs_early = ifelse(Study.Group %in% c("HS", "Stage_I_II", "Stage_III_IV"), NA, 
                            ifelse(Study.Group %in% c("Stage_0", "MP"), "Stage_0_MP", "Healthy"))
  )

#add m__ before each metabolite. 
#this will allow us to parse and automate the data downstream
metabolites <- metabolites %>%
  rename_with(~ paste0("m__", .), -Sample)
# Apply log2 transformation only to numeric columns, excluding "Sample"
metabolites <- metabolites %>%
  mutate(across(where(is.numeric) & !all_of("Sample"), ~ log2(. + 1)))

#merge species and genus by Sample
mss <- merge(genera, species, by = "Sample")
# Calculate sums of numeric columns excluding 'Sample'
sums <- mss %>%
  summarise(across(where(is.numeric) & !all_of("Sample"), sum, na.rm = TRUE))
# Select columns with sums >= 0.001 excluding 'Sample'
filtered_mss <- mss %>%
  select(Sample, where(~ !is.numeric(.) || sum(., na.rm = TRUE) >= 0.001))
# Calculate sums again on the filtered data, excluding 'Sample'
filtered_sums <- filtered_mss %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


#add t__ before all taxa to make easier to parse
filtered_mss<- filtered_mss %>%
  rename_with(~ paste0("t__", .), -Sample)


summary(metadata$`Brinkman Index`)
summary(metadata$Alcohol)

#log transform Fecal.calprotectin
metadata$`Brinkman Index` <- log2(metadata$`Brinkman Index`+1)
metadata$Alcohol <- log2(metadata$Alcohol+1)


metadata_to_merge <- metadata[, ..metadata_cols]

#check NAs
table(metadata_to_merge$Gender)
summary(metadata_to_merge$Age)
summary(metadata_to_merge$BMI)
summary(metadata_to_merge$`Brinkman Index`)
summary(metadata_to_merge$Alcohol)

#merge metadata and metabolites
metabolites_with_labels <- merge(metadata_to_merge, metabolites, by = "Sample")
yachida_processed_data <- merge(metabolites_with_labels, filtered_mss, by = "Sample")


##THIS SECTION IS TO GENERATE non-processed datasets (no differential expression - as controls).
yachida_data<- yachida_processed_data
yachida_gender_age_bmi_brinkman_alcohol <- yachida_data
yachida_healthyvscancer <- yachida_data  %>% filter(!is.na(disease_status))
yachida_healthyvsearly <- yachida_data  %>% filter(!is.na(healthy_vs_early))
yachida_healthyvsstageI_II <- yachida_data %>% filter(Study.Group %in% c("Healthy", "Stage_I_II"))
yachida_healthyvsstageIII_IV <- yachida_data %>% filter(Study.Group %in% c("Healthy", "Stage_III_IV"))

table(yachida_data$Study.Group)

write.csv(yachida_gender_age_bmi_brinkman_alcohol, file = "nofilter_yachida_gender_age_bmi_brinkman_alcohol_12062024.csv", row.names = FALSE) 
write.csv(yachida_healthyvscancer, file = "nofilter_yachida_healthyvscancer_12062024.csv", row.names = FALSE) 
write.csv(yachida_healthyvsearly, file = "nofilter_yachida_healthyvsearly_12062024.csv", row.names = FALSE) 
write.csv(yachida_healthyvsstageI_II, file = "nofilter_yachida_healthyvsstageI_II_12062024.csv", row.names = FALSE) 
write.csv(yachida_healthyvsstageIII_IV, file = "nofilter_yachida_healthyvsstageIII_IV_12062024.csv", row.names = FALSE) 


#SECTION 2: DIFFERENTIAL EXPRESSION ANALYSIS FOR FEATURE REDUCTION. 
#will perform separately for each dataset. 

#PERFORM LIMMA FOR METABOLITES
m_columns <- grep("^m__", colnames(yachida_processed_data), value = TRUE)
metabolite_data <- yachida_processed_data %>% select(m_columns)
rownames(metabolite_data) <- yachida_processed_data$Sample
metabolite_data <- as.matrix(metabolite_data)
metabolite_data_t <- t(metabolite_data)
metadata_metab <- yachida_processed_data %>% 
  select(Study.Group) %>% 
  as.data.frame()
# Set row names using the Sample column from yachida_processed_data
rownames(metadata_metab) <- yachida_processed_data$Sample

#use metadata, design and contrasts from above for this problem. 
#contrasts fit
# Fit the linear model
design <- model.matrix(~ Study.Group-1, data = metadata_metab)
# Define contrasts
contrasts <- makeContrasts(
  HealthyvsMP = Study.GroupHealthy - Study.GroupMP,
  HealthyvsStage_0 = Study.GroupHealthy - Study.GroupStage_0,
  HealthyvsStage_I_II = Study.GroupHealthy - Study.GroupStage_I_II,
  HealthyvsStage_III_IV = Study.GroupHealthy - Study.GroupStage_III_IV,
  levels = design
)
fit_m <- lmFit(metabolite_data_t, design)
# Apply contrasts
fit_contrasts_m <- contrasts.fit(fit_m, contrasts)
# Apply empirical Bayes moderation
fit_contrasts_m <- eBayes(fit_contrasts_m)


# Define a function to extract and filter results based on adjusted p-value and fold change
extract_significant_taxa <- function(fit_contrasts, contrast_name, p_value_threshold = 0.05, fold_change_threshold = 2.0, filter_results = FALSE) {
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
  "HealthyvsMP",
  "HealthyvsStage_0",
  "HealthyvsStage_I_II",
  "HealthyvsStage_III_IV"
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




# # # Define a function to extract and filter results based on adjusted p-value and fold change
# extract_significant_taxa <- function(fit_contrasts, contrast_name, p_value_threshold = 0.05, fold_change_threshold = 2.0) {
#   # Extract results for the given contrast
#   results <- topTable(fit_contrasts, coef = contrast_name, adjust = "BH", p.value = p_value_threshold, number = Inf)
#   # Filter results based on adjusted p-value and fold change
#   significant_hits <- results[
#     results$adj.P.Val < p_value_threshold &
#       (results$logFC > fold_change_threshold | results$logFC < -fold_change_threshold),
#   ]
#   return(significant_hits)
# }



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
saveWorkbook(wb, file = "Yachida_Differential_Abundance_Results_with_Metabolites_2.0fc_0.05pval_12062024.xlsx", overwrite = TRUE)

#these will be uncoupled in three ways based on ControlvDisease, ControlvUC and ControlvCD. 
#reduce features to significant only. We are going to use these to generate some models.
significant_metabolites <- lapply(results_list_m, function(df) {
  rownames_list <- rownames(df)
  return(rownames(df))
})
#this will be used for the metabolomics data. 
diff_expressed_metabolites_all <- unique(unlist(significant_metabolites, recursive = FALSE))
#only two in total. Tried 1.5 and was 4. Will just use all the metabolites for this analysis. 


#PERFORM WILCOXON FOR TAXA
t_columns <- grep("^t__", colnames(yachida_processed_data), value = TRUE)
w_mss <- yachida_processed_data %>% select('Study.Group', t_columns)
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
control_group <- "Healthy"
treatment_groups <- c("MP", "Stage_0", "Stage_I_II", "Stage_III_IV")

table(metadata$Study.Group)

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
#results <- results[!is.na(results$adjusted_p) & results$adjusted_p <= 0.05, ]

split_data <- split(results, results$Comparison)
names(split_data) <- gsub("/", "_", names(split_data))
# Specify the output Excel file path
output_file <- "yachida_Differential_Abundance_Results_with_Taxa_0.05pval_12062024.xlsx"
# Write each dataframe to a separate sheet in the Excel file
write_xlsx(split_data, path = output_file)

#read in the wilcoxon data
library(readxl)  
file_name <- "yachida_Differential_Abundance_Results_with_Taxa_0.05pval_12062024.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(file_name)

# Read all sheets into a list
mss_results <- lapply(sheet_names, function(sheet) {
  read_excel(file_name, sheet = sheet)
})
# Assign sheet names to the list elements
names(mss_results) <- sheet_names
#results <- results[!is.na(results$p_value) & results$p_value <= 0.05, ]

#keep only unique values for df1 and df2
df1 <- mss_results[[1]]
df2 <- mss_results[[2]]
# Combine column names from both dataframes and keep only unique values
unique_columns_MPStage0 <- unique(c(df1$Taxon, df2$Taxon))

df3 <- mss_results[[3]]
df4 <- mss_results[[4]]
sum(df4$adjusted_p < 0.05, na.rm = TRUE)

#get unique taxa for all_df
all_unique_taxa <- unique(unlist(lapply(mss_results, function(df) {
  df$Taxon
})))



#SECTION 3: COMBINE ABX METAB AND TAXA - COMBINE DIET METAB AND TAXA

#we will use all metabolite data since little diff expressed data even with p-vals. 
diff_expressed_data_list_all <- c(colnames(metabolites_with_labels), all_unique_taxa)
yachida_all_data <- yachida_processed_data[, ..diff_expressed_data_list_all]
# yachida_healthyvsMP <- c(colnames(metabolites_with_labels), colnames(df1))
# yachida_healthyvsMP_data <- yachida_processed_data[, ..yachida_healthyvsMP]
yachida_healthyvsstage0_MP <- c(colnames(metabolites_with_labels), unique_columns_MPStage0)
yachida_healthyvsstage0_MP_data <- yachida_processed_data[, ..yachida_healthyvsstage0_MP]
yachida_healthyvsstage_I_II <- c(colnames(metabolites_with_labels), df3$Taxon)
yachida_healthyvsstageI_II_data <- yachida_processed_data[, ..yachida_healthyvsstage_I_II]
yachida_healthyvsstage_III_IV <- c(colnames(metabolites_with_labels), df4$Taxon)
yachida_healthyvsstageIII_IV_data <- yachida_processed_data[, ..yachida_healthyvsstage_III_IV]

yachida_healthyvscancer_filter <- yachida_all_data  %>% filter(!is.na(disease_status))
yachida_healthyvsearly_filter <- yachida_healthyvsstage0_MP_data  %>% filter(!is.na(healthy_vs_early))
yachida_healthyvsstageI_II_filter <- yachida_healthyvsstageI_II_data %>% filter(Study.Group %in% c("Healthy", "Stage_I_II"))
yachida_healthyvsstageIII_IV_filter <- yachida_healthyvsstageIII_IV_data %>% filter(Study.Group %in% c("Healthy", "Stage_III_IV"))

write.csv(yachida_all_data, file = "filter_yachida_gender_age_bmi_brinkman_alcohol_12092024.csv", row.names = FALSE) 
write.csv(yachida_healthyvscancer_filter, file = "filter_yachida_healthyvscancer_12092024.csv", row.names = FALSE) 
write.csv(yachida_healthyvsearly_filter, file = "filter_yachida_healthyvsearly_12092024.csv", row.names = FALSE) 
write.csv(yachida_healthyvsstageI_II_filter, file = "filter_yachida_healthyvsstageI_II_12092024.csv", row.names = FALSE) 
write.csv(yachida_healthyvsstageIII_IV_filter, file = "filter_yachida_healthyvsstageIII_IV_12092024.csv", row.names = FALSE) 