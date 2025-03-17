#!/bin/bash

# Base directory containing the configuration files
#change lines below to run bash script
base_dir="/home2/s180020/Desktop/microbiome-metabolome-curated-data/IntegratedLearner"
config_dir="$base_dir/ReRun_Configs"

# Path to the job template script
job_template="$base_dir/job_template.sh"

# Loop through each subdirectory in the Configuration_Files directory
for subdir in "$config_dir"/*/; do
  # Remove trailing slash
  subdir=${subdir%/}
  
  #echo "Processing subdirectory: $subdir"

  # Loop through each configuration file in the subdirectory
  for config_file in "$subdir"/*.txt; do
    # Extract just the file name without the path
    config_file_name=$(basename "$config_file")
    
    #echo "Submitting job for $config_file_name in subdirectory $subdir"

    # Submit the sbatch job with the configuration file as an argument
    source $config_file
    study_name=${study_name//\'/}

#    study_name=`cat $config_file |grep study_name |cut -f2 -d'='`
    if [ -e "results/AnalysisRunDocumentation/${study_name}.txt" ]; then
      if grep -q Elapsed "results/AnalysisRunDocumentation/${study_name}.txt"  ; then
        echo "yes,$study_name";
      else
       echo "no,$study_name";
      fi
    else
       echo "no_result,$study_name"
    fi
  done
done
