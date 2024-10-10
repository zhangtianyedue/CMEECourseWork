#!/bin/bash
# Author: zhangtianye
# Script: convert_images.sh
# Description: Converts all files with a specific extension (e.g., .tif) to another format (e.g., .png).
# Arguments: 1 -> Input file extension (e.g., tif)
#            2 -> Output file extension (e.g., png)
# Date: Oct 2024

#!/bin/bash

# Check if both input and output directories are provided
if [ $# -eq 2 ]; then
    input_dir=$1
    output_dir=$2

    # Check if input directory exists
    if [ ! -d "$input_dir" ]; then
        echo "Error: Input directory $input_dir does not exist."
        exit 1
    fi

    # Check if output directory exists, if not, create it
    if [ ! -d "$output_dir" ]; then
        echo "Output directory $output_dir does not exist. Creating it..."
        mkdir -p "$output_dir"
    fi

    # Loop through all .tif files in the input directory and convert them
    for f in "$input_dir"/*.tif; do
        if [ -f "$f" ]; then
            echo "Converting $f"
            output_file="$output_dir/$(basename "$f" .tif).png"
            convert "$f" "$output_file"
        else
            echo "No .tif files found in $input_dir."
        fi
    done

# If no input or output directories are provided, execute the original code
else
    echo "No input or output directories specified. Running in the current directory..."

    # Original code to convert .tif files in the current directory
    for f in *.tif; do
        echo "Converting $f"
        convert "$f" "$(basename "$f" .tif).png"
    done
fi
