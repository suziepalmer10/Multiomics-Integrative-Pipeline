# Multiomics-Integrative-Pipeline
=======
# Integrative Pipeline Code Instructions


## Table of Contents
- [SECTION 1: Manuscript Code](#section-1-manuscript-code)  
  - [1) Original Studies Used in the Analysis](#1-original-studies-used-in-the-analysis)  
  - [2) Multi-omics Data Processing by Borenstein Lab](#2-multi-omics-data-processing-by-borenstein-lab)  
  - [3) Pre-processing for the Integrative Pipeline](#3-pre-processing-for-the-integrative-pipeline)  
  - [4) Script Storage for Integrative Analysis](#4-script-storage-for-integrative-analysis)  
  - [5) Running the Integrative Analysis](#5-running-the-integrative-analysis)  
  - [6) Processing the Performance Data for the Tables Used for Heatmap](#6-processing-the-performance-data-for-the-tables-used-for-heatmap)  
  - [7) Combined Heatmap Creation](#7-combined-heatmap-creation)  
  - [8) Pie Charts](#8-pie-charts)  
  - [9) Feature Extraction for Model Comparisons](#9-feature-extraction-for-model-comparisons)  
  - [10) Elastic Net, Random Forest, and XGBoost Feature Comparison Script and Violin Plot Creation](#10-elastic-net-random-forest-and-xgboost-feature-comparison-script-and-violin-plot-creation)  
  - [Citation](#citation)  
  - [License](#license)  

---

## SECTION 1: Manuscript Code

### 1) Original Studies Used in the Analysis
The following studies provided multi-omics microbiome data:

- **Yachida_CRC_2019**: Yachida, Shinichi, et al. "Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the gut microbiota in colorectal cancer." *Nature Medicine* 25.6 (2019): 968-976.
- **Erawijantari_Gastrectomy_2020**: Erawijantari et al. "Influence of gastrectomy for gastric cancer treatment on faecal microbiome and metabolome profiles." *Gut* 69.8 (2020): 1404-1415.
- **Wang_ESRD_2020**: Wang, Xifan, et al. "Aberrant gut microbiota alters host metabolome and impacts renal failure in humans and rodents." *Gut* 69.12 (2020): 2131-2142.
- **Franzosa_IBD_2019**: Franzosa, Eric A., et al. "Gut microbiome structure and metabolic activity in inflammatory bowel disease." *Nature Microbiology* 4.2 (2019): 293-305.

### 2) Multi-omics Data Processing by Borenstein Lab

- The datasets from Yachida, Erawijantari, Wang, and Franzosa were obtained from the [Borenstein Lab repository](https://github.com/borenstein-lab/microbiome-metabolome-curated-data/wiki/Data-overview#datasets-included).
- Notes on metagenomic shotgun sequencing data processing can be found [here](https://github.com/borenstein-lab/microbiome-metabolome-curated-data/wiki/Data-processing-details#metagenomics-processing-notes). Taxonomic classification was performed using Kraken v2.1.1 and Bracken v2.8.
- Notes on metabolomics data processing can be found [here](https://github.com/borenstein-lab/microbiome-metabolome-curated-data/wiki/Data-processing-details#metabolomics-processing-notes).

### 3) Pre-processing for the Integrative Pipeline

> ⚠️ **Important:** Ensure correct file paths before running the pipeline to avoid errors.

#### Data Tables Used
- **Metagenomics**: `genera.tsv`, `species.tsv`
- **Metabolomics**: `mtb.tsv`
- **Metadata**: `metadata.tsv`

#### Data Processing Steps
**Metagenomics:**
- `genera.tsv` and `species.tsv` were combined.
- "t__" was prefixed to taxa names for easier parsing.
- Taxa with abundance <0.001 across all samples were removed.
- Wilcoxon test was performed for differential abundance analysis (*p*-value < 0.05).

**Metabolomics:**
- "m__" was prefixed to metabolite names for easier parsing.
- Log scaling was applied based on dataset distribution.
- Limma analysis was performed with log-fold-change >2 and/or *p*-value < 0.05 as thresholds.

**Metadata:**
- Missing values for response variables were removed.
- Log scaling was applied to response variables (e.g., BMI) based on distribution.

**Data Merging:**
- Metagenomics, metabolomics, and metadata were merged based on the `Sample` metadata column.

**Scripts for Data Processing**
- Scripts for this code can be found in: `Manuscript_Code/A_Data_Processing/`. There are 4 directories contained inside: Erawijantari, Wang, Franzosa and Yachida. 
- Each of these directories contain two subdirectories: `/A_Data_ObtainedfromBorenstein` and `/B_OutputDataForIntegrativeAnalysis`. 
- The scripts to process each dataset are also contained in these directories: `*_data_processing.R` (i.e., `Erawijantari_data_processing.R`) and are used to uniquely process each study.  
- `/A_Data_ObtainedfromBorenstein`: contains the originally processed files from Borenstein: `species.tsv`, `genera.tsv`, `mtb.tsv` and `metadata.tsv`. 
- `/B_OutputDataForIntegrativeAnalysis`: Contains the output files from running the `*_data_processing.R` file. Notice that there will be more than one output file based on missing values in response variables (which can change the sample size) and whether the data is reduced or full-dimensional. 

### 4) Script Storage for Integrative Analysis
- The scripts to run the Integrative Pipeline can be found at `Manuscript_Code/B_Integrative_Pipeline/A_Integrative_Pipeline_Code`. There are three subdirectories:
- `A_Integrative_Pipeline_Scripts`: this contains the scripts needed to run the Integrative Pipeline. 
  * `automated_pipeline.R` is used as the main script and feeds in data, desired model (elastic net (`model_functions/enet_function.R`), random forest (`model_functions/rf_function.R`), or xgboost (`model_functions/xgboost_function.R`) and either binary function files (`model_functions/binary_functions.R` and `model_functions/binary_performance_functions.R`) or continuous functions (`model_functions/continuous_functions.R` and `model_functions/continuous_performance_functions.R`). 
  * To modify the grids for any of the machine learning models above, you must go to the separate script (i.e., `model_functions/enet.R`) . 
  * Performance metrics for continuous are MAE, RMSE and R2.
  * Performance metrics for binary are Accuracy, Kappa and AUROC.
- `B_BashScripts`: `job_template.sh`, `submit_jobs.sh` and `submit_jobs_individually.sh` are the three scripts associated with running all data associated with the metadata files. The `job_template.sh` will feed the metadata.txt script (shown above) and associated parameters into the `automated_pipeline.R` script and ensures that the correct performance metrics and visualization scripts are read in. The `submit_jobs.sh` script will iterate through all directories contained in the configuration directory (config_dir) and submit all metadata files to the cluster. The `submit_jobs_individually.sh` will only submit one directory of metadata files to run on the cluster. `check_jobs.sh` is used to quickly check whether jobs have finished or not by using the `Analysis_Run_Documentation` file. 
- `C_Configuration_Files`: contains the configuration files for all of the analyses used for this study. 

### 5) Running the Integrative Analysis
- You will need to fill out configuration files (or use the ones contained in `C_Configuration_Files`. 
- The parameters need for each configuration file are:
```bash
model_to_run= <char> 
file_path= <char> 
input_file= <char> 
study_name= <char> 
type_of_analysis= "continuous" or "binary"
response_variable= <char> 
stratify_variable= <char> 
training_proportion= <int> 
num_repeats= <int>
num_folds= <int> 
```
> ⚠️ **Important:** Note: if any of these parameters are not filled out, the `automated_pipeline.R` script will throw and error. You can leave 'stratfy_variable' blank with '' if you do not wanted the data stratified. Make sure all paths are correct. 

The `job_template.sh` is used to parse the information provided from the configuration files. 
The `B_BashScripts/submit_jobs_individual.sh` will submit an individual directory. 
The `B_BashScripts/submit_jobs.sh` will recursively submit all subdirectory files contained in a directory. 
Once the user has the correct configuration files and correct paths, they can rerun these analyses using either of the `submit_jobs_*.sh` scripts. 

```bash
bash B_BashScripts/submit_jobs.sh
```
Each dataset run generates the following output files:

| Output File Type           | Description                                         |
|----------------------------|-----------------------------------------------------|
| `performance_csv`          | Performance metrics in table format                |
| `Analysis_Run_Documentation` | Run details, file information, etc.                  |
| `Visualizations`           | PNG images of performance metrics                  |
| `RDataFiles`               | Saved RData image for each run (*not stored on GitHub*) |

Performance CSV files are stored in `Manuscript_Code/B_Integrative_Pipeline/B_Output_PerformanceMetrics` by study.

### 6) Processing the Performance Data for the Tables Used for Heatmap

Under `Manuscript_Code/B_Integrative_Pipeline/C_Output_AUROC_RMSE_Heatmap_Tables` are the performance output files contained in subdirectories for each dataset (i.e., Erawijantari).
  * Binary and Continuous validation and performance metrics for feature reduced and full dimensional datasets are contained in each subdirectory (8 in total). 
  * The `*_performance_heatmaps.R` (i.e., `Erawijantari_performance_heatmaps.R`) are contained for each dataset for extracting and reformatting the data. 

### 7) Combined Heatmap Creation

- Under `Manuscript_Code/B_Integrative_Pipeline/D_Heatmap`, the code for combining and creating the heatmaps for Figure 2-3 and Supplemental Figure 1 and 2. This can be replicated using the `Combinedheatmaps_12152024.R`.
- The Heatmap images are provided as pdf files in this directory. 

### 8) Pie Charts

- Pie Charts that appear Figure 2-3 and Supplemental Figure 1 and 2 are contained under `Manuscript_Code/B_Integrative_Pipeline/E_PieCharts`. These represent the top 3 performing models for each model (Elastic Net, Random Forest and XGBoost) created for a response variable. 
- The PieChart images are provided as pdf files in this directory. 

## 9) Feature Extraction for Model Comparisons

- Files for this analysis are contained in: `Manuscript_Code/C_Feature_Importance_Values`. 
- The R Markdown file: `Feature_extraction.Rmd` or Rscript: `feature_extraction_12302024.R` can be used to extract features. 
- This analysis extracts features from Elastic Net, Random Forest and XGBoost models. The absolute value of the features was taken and Min-Max Normalization was applied to enable comparability between the model feature importances. 
- The Output for this Feature Extraction Analysis is contained in a subdirectory: `Output_FeatureExtraction`. 

## 10) Elastic Net, Random Forest and XGBoost Feature Comparison Script and Violin Plot Creation

- Files for this analysis are contained in: `Manuscript_Code/D_Feature_Comparisons`.
- The R Markdown file: `FeatureExtractionAnalysis.Rmd` is used to compare different proportions of features between the Elastic Net, Random Forest and XGBoost models, which are used for Figures 4-6. 
- The default proportion compared between the models are: 1, 5, 10, 20 and 50 and this is stored in the subdirectory: `Output_FeatureComparison`. 
- The R Markdown file: `Feature_Comparison_ViolinPlot.Rmd` is used to create Violin Plots to visualize model feature selection comparisons for Elastic Net, Random Forest and XGBoost models, which are used for Figures 4-6. 
- The output Violin plots are stored in the subdirectory: `Output_ViolinPlots`.

---

## 🔬 SECTION 2: Running the Integrated Learner Pipeline on Your Own Data

If you are interested in running the Integrated Learner Pipeline on your own dataset, navigate to `Code_and_Analyses/RunYourOwnData`. Follow the same setup and configuration procedures outlined in **SECTION 1**, ensuring your dataset adheres to the required format.

For additional assistance, please refer to the documentation in `RunYourOwnData`.

---

## 📖 Citation
If you use this pipeline in your research, please cite: 

---

## 📜 License
This project is open-source under the MIT License.

---

For any issues or questions, please open an issue on GitHub and/or contact suzette.palmer@utsouthwestern.edu.


