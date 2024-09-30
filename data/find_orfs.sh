#!/bin/bash

set -x

# Ensure the script takes two arguments
input_directory="/home/mullara/GeneFinder/GeneFinderTool/data"
output_directory="$1"  # Output directory is passed as the first argument to this script

# Check if the output directory argument is supplied
if [ -z "$output_directory" ]; then
    echo "Usage: $0 <output_directory>"
    exit 1
fi

# Ensure the output directory exists
mkdir -p "$output_directory"

find "$input_directory" -type f -name "*.fna" | while read -r file; do
    echo "Processing $file..."
    python3 process_orfs.py "$file" "$output_directory"
done
