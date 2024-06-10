#!/bin/bash
# Description: this script uses Qualimap RNASEQ and RSeQC read_distribution.py to assess  and  distribution of mapped reads
# Usage: ./mapping_check.sh

# input file includes name of samples
samples="samples.txt"
if [ ! -f "$samples" ]; then
    echo "Input file ${samples} not found."
    exit 1
fi

# create output directory to save output
mkdir -p ./mapping_check/qualimap/
mkdir -p ./mapping_check/rseqc/

# define log file path
log_file_path_qualimap="./mapping_check/qualimap/rnaseq.log"
log_file_path_rseqc="./mapping_check/rseqc/rnaseq.log"

# start logging QUALIMAP RNASEQ
echo "$(date) - Starting QUALIMAP RNASEQ" 2>&1 | tee -a $log_file_path_qualimap
qualimap --version 2>&1 | tee -a $log_file_path_qualimap
printf "\n\n\n" 2>&1 | tee -a $log_file_path_qualimap

# QUALIMAP RNASEQ: loop through each sample
while IFS= read -r sample; do
    echo "$(date) - Running qualimap rnaseq for ${sample}" 2>&1 | tee -a $log_file_path_qualimap
    qualimap rnaseq -bam ./processed/${sample}_sorted.bam -gtf ./ref/gh38.gtf -outformat HTML -outfile ${sample} -outdir ./mapping_check/qualimap/${sample} --paired --sorted --java-mem-size=4G -npb 50 2>&1 | tee -a $log_file_path_qualimap
    echo "$(date) - Completed converting for ${sample}" 2>&1 | tee -a $log_file_path_qualimap
    printf "\n" 2>&1 | tee -a $log_file_path_qualimap

done < "$samples"

echo "$(date) - Completed QUALIMAP RNASEQ" 2>&1 | tee -a $log_file_path_qualimap

# merge output
echo "$(date) - Merging QUALIMAP RNASEQ output" 2>&1 | tee -a $log_file_path_qualimap
multiqc ./mapping_check/qualimap/
echo "$(date) - Merging completed" 2>&1 | tee -a $log_file_path_qualimap



# start logging RSeQC
echo "$(date) - Starting RSeQC read_distribution.py" 2>&1 | tee -a $log_file_path_rseqc

# RSeQC read_distribution.py: loop through each sample
while IFS= read -r sample; do
    out_f="./mapping_check/rseqc/${sample}.txt"
    echo "$(date) - Assessing read distribution for ${sample}" 2>&1 | tee -a $log_file_path_rseqc
    read_distribution.py -i ./mapped/${sample}_sorted.bam -r ./ref/gh38.bed12 2>&1 | tee -a $log_file_path_rseqc
    echo "$(date) - Completed assessing ${sample}" 2>&1 | tee -a $log_file_path_rseqc | tee -a $out_f
    printf '\n' 2>&1 | tee -a $log_file_path_rseqc

done < "$samples"

echo "$(date) - Completed RSeQC" 2>&1 | tee -a $log_file_path_rseqc

# merge output
echo "$(date) - Merging RSeQC read_distribution.py output" 2>&1 | tee -a $log_file_path_rseqc
multiqc ./mapping_check/rseqc/*.txt
echo "$(date) - Merging completed" 2>&1 | tee -a $log_file_path_rseqc
