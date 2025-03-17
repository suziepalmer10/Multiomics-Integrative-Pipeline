#!/bin/bash

cd /home2/s180020/Desktop/microbiome-metabolome-curated-data/IntegratedLearner
# Directory containing the configuration files
#config_dir="Configuration_Files/Redo-jobscancelled"
config_dir="TryAgainSucker"


# Loop through each configuration file in the directory
for config_file in "$config_dir"/*.txt; do
  echo "Submitting job for $config_file"
  
  # Submit the sbatch job with the configuration file as an argument
  sbatch -p 128GB job_template.sh "$config_file"

done
