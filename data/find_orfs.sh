#!/bin/bash

set -x

find /home/mullara/GeneFinder/GeneFinderTool/data -type f -name "*.fna" | while read -r file; do
    echo "Processing $file..."
    python3 process_orfs.py "$file"
done


