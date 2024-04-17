#!/bin/bash
# Description: this script checks the quality of reads in FASTQ files
# Usage: ./fastqc.sh

# input file includes name of samples
samples="samples.txt"
if [ ! -f "$samples" ]; then
    echo "Input file '$samples' not found."
    exit 1
fi

# create output directory to save quality results
mkdir -p ./fastqc

# define log file path
log_file_path="./fastqc/fastqc.log"

# start logging
echo "$(date) - Running FASTQC" 2>&1 | tee -a $log_file_path
fastqc -v 2>&1 | tee -a $log_file_path
printf '\n\n\n' 2>&1 | tee -a "$log_file_path"

# loop through each sample
while IFS= read -r sample; do

    echo "$(date) - Processing ${sample}" 2>&1 | tee -a $log_file_path
    fastqc "./fastq/${sample}_1.fq.gz" "./fastq/${sample}_2.fq.gz" -o ./fastqc
    echo "$(date) - Completed processing ${sample}" 2>&1 | tee -a $log_file_path
    printf '\n' 2>&1 | tee -a "$log_file_path"

done < "$samples"

echo "$(date) - Completed FASTQC" 2>&1 | tee -a $log_file_path

echo "$(date) - Merging results using MULTIQC" 2>&1 | tee -a $log_file_path
multiqc ./fastqc/ 2>&1 | tee -a $log_file_path
echo "$(date) - Merging completed" 2>&1 | tee -a $log_file_path

