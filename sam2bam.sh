#!/bin/bash
# Description: this script uses samtools to convert SAM (from HISAT2) to BAM, sort, coordinate, index,...
# Usage: ./samtools.sh

# input file includes name of samples
samples="samples.txt"
if [ ! -f "$samples" ]; then
    echo "Input file '$samples' not found."
    exit 1
fi

# create output directory to save .bam files
mkdir -p ./processed

# define log file path
log_file_path="./processed/samtools.log"

# start logging
echo "$(date) - Running Samtools" 2>&1 | tee -a $log_file_path
samtools --version 2>&1 | tee -a $log_file_path
printf '\n\n\n' 2>&1 | tee -a "$log_file_path"

# loop through each sample
while IFS= read -r sample; do

    # SAM to BAM
    echo "$(date) - Converting SAM file: ${sample}" 2>&1 | tee -a $log_file_path
    samtools view -bS "./mapped/${sample}.sam" > "./processed/${sample}.bam" 2>&1 | tee -a $log_file_path
    echo "$(date) - Completed converting for ${sample}" 2>&1 | tee -a $log_file_path
    printf '\n' 2>&1 | tee -a "$log_file_path"
    
    # Sorting BAM
    echo "$(date) - Sorting BAM file: ${sample}" 2>&1 | tee -a $log_file_path
    samtools sort -@ 3 "./mapped/${sample}.bam" -o "./mapped/${sample}_sorted.bam" 2>&1 | tee -a $log_file_path
    echo "$(date) - Sorting completed for ${sample}" 2>&1 | tee -a $log_file_path
    printf '\n' 2>&1 | tee -a "$log_file_path"

    # Indexing BAM
    echo "$(date) - Indexing BAM file: ${sample}" 2>&1 | tee -a $log_file_path
    samtools index -@ 3 "./mapped/${sample}_sorted.bam" 2>&1 | tee -a $log_file_path 2>&1 | tee -a $log_file_path
    echo "$(date) - Indexing completed for ${sample}" 2>&1 | tee -a $log_file_path
    printf '\n' 2>&1 | tee -a "$log_file_path"

done < "$samples"

echo "$(date) - Completed SAMTOOLS" 2>&1 | tee -a $log_file_path
