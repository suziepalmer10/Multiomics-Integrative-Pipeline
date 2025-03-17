#!/bin/bash

# Source the configuration file
echo $1
cat $1
source $1

# Remove surrounding quotes from the values
model_to_run=${model_to_run//\'/}
file_path=${file_path//\'/}
input_file=${input_file//\'/}
study_name=${study_name//\'/}
type_of_analysis=${type_of_analysis//\'/}
response_variable=${response_variable//\'/}
stratify_variable=${stratify_variable//\'/}
training_proportion=${training_proportion//\'/}
num_repeats=${num_repeats//\'/}
num_folds=${num_folds//\'/}

#SBATCH --job-name=$study_name     # Job name
#SBATCH --output=${study_name}.log
#SBATCH --error=${study_name}.err
#SBATCH --partition=32GB

# Send an email when the job status changes, to the specified address.
#SBATCH --mail-type ALL
#SBATCH --mail-user suzette.palmer@utsouthwestern.edu

module load python/3.8.x-anaconda
conda activate R4.2
cd /home2/s180020/Desktop/microbiome-metabolome-curated-data/IntegratedLearner

#defined arguments to use
Rscript automated_pipeline.R \
--model_to_run "$model_to_run" \
--file_path "$file_path" \
--input_file "$input_file" \
--study_name "$study_name" \
--type_of_analysis "$type_of_analysis" \
--response_variable "$response_variable" \
--stratify_variable "$stratify_variable" \
--training_proportion "$training_proportion" \
--num_repeats "$num_repeats" --num_folds "$num_folds"
