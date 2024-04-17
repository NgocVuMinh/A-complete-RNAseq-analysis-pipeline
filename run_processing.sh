#!/bin/bash
# Usage: ./run_processing.sh

# ./download.sh
./quality_control.sh
./mapping.sh
./sam2bam.sh
./mapping_check.sh
./quantification.sh