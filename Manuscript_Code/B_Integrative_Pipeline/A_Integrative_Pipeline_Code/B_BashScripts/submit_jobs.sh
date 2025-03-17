#!/bin/bash

# Base directory containing the configuration files
base_dir="/home2/s180020/Desktop/microbiome-metabolome-curated-data/IntegratedLearner"
#change line below for directory
config_dir="$base_dir/yachida_ConfigFiles"

# Path to the job template script
job_template="$base_dir/job_template.sh"

# Loop through each subdirectory in the Configuration_Files directory
for subdir in "$config_dir"/*/; do
  # Remove trailing slash
  subdir=${subdir%/}
  
  echo "Processing subdirectory: $subdir"

  # Loop through each configuration file in the subdirectory
  for config_file in "$subdir"/*.txt; do
    # Extract just the file name without the path
    config_file_name=$(basename "$config_file")
    
    echo "Submitting job for $config_file_name in subdirectory $subdir"

    # Submit the sbatch job with the configuration file as an argument
    sbatch -p 32GB --job-name="$(basename "$subdir")_$(basename "$config_file" .txt)" "$job_template" "$config_file"
  done
done
