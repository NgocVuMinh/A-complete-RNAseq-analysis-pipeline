#!/bin/bash
# Description: this script uses featureCounts to quantify genes using BAM files
# Usage: ./featureCounts.sh

# input file includes name of samples
samples="samples.txt"
if [ ! -f "$samples" ]; then
    echo "Input file '$samples' not found."
    exit 1
fi

# create output directory to save gene count results
mkdir -p ./featureCounts

# define log file path
log_file_path="./featureCounts/featureCounts.log"

# start logging
echo "$(date) - Running featureCounts" 2>&1 | tee -a $log_file_path
featureCounts -v 2>&1 | tee -a $log_file_path
printf '\n\n\n' 2>&1 | tee -a "$log_file_path"

annotation="./ref/gh38.gtf"

# loop through each sample
while IFS= read -r sample; do
  
    echo "$(date) - Quantifying genes for ${sample}" 2>&1 | tee -a $log_file_path
    featureCounts -p -s 0 -T 3 -a $annotation -o "./featureCounts/${sample}_featureCounts.txt" "./processed/${sample}.bam" --verbose 2>&1 | grep -v "is marked as paired, but its mate does not occur next to it in your BAM file" | tee -a $log_file_path 
    cut -f1,7,8,9,10,11,12 ./featureCounts/${sample}_featureCounts.txt > ./featureCounts/${sample}_featureCounts_clean.txt
    echo "$(date) - Completed quantifying genes for ${sample}" 2>&1 | tee -a $log_file_path
    printf '\n' 2>&1 | tee -a "$log_file_path"

done < "$samples"

echo "$(date) - Completed featureCounts" 2>&1 | tee -a $log_file_path