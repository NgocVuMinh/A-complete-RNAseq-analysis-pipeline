#!/bin/bash
# Description: this script uses HISAT2 to map reads to reference genome
# Usage: ./hisat2.sh

# input file includes name of samples
samples="samples.txt"
if [ ! -f "$samples" ]; then
    echo "Input file '$samples' not found."
    exit 1
fi

# create output directory to save .sam files
mkdir -p ./mapped
mkdir -p ./mapped/report
 
# define log file path
log_file_path="./mapped/report/hisat2_mapping.log"

# start logging
echo "$(date) - Starting mapping with HISAT2" 2>&1 | tee -a $log_file_path
hisat2 --version 2>&1 | tee -a $log_file_path
printf "\n\n\n" 2>&1 | tee -a "$log_file_path"

# loop through each sample
while IFS= read -r sample; do
    echo "$(date) - Mapping ${sample}" 2>&1 | tee -a $log_file_path
    hisat2 -x "./ref/gh38_index/genome" -1 "./fastq/${sample}_1.fq.gz" -2 "./fastq/${sample}_2.fq.gz" -S "./mapped/${sample}.sam" --summary-file "./mapped/report/${sample}_report.txt" --time --threads 3 2>&1 | tee -a $log_file_path
    echo "$(date) - Completed mapping ${sample}" 2>&1 | tee -a $log_file_path
    printf "\n" 2>&1 | tee -a "$log_file_path"

done < "$samples"

echo "$(date) - Completed HISAT2" 2>&1 | tee -a $log_file_path
