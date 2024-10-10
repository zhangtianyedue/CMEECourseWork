#!/bin/bash
# Author: zhangtianye
# Script: countlines.sh
# Description: Counts the number of lines in one or more files and prints the result.
# Arguments: 1 or more file names
# Date: Oct 2024

# Check if at least one argument (file name) is provided
if [ $# -lt 1 ]; then
    echo "Error: Please provide at least one file name."
    echo "Usage: ./countlines.sh <file1> [file2 ...]"
    exit 1
fi

# Loop through all provided files
for file in "$@"; do
    # Check if the file exists
    if [ ! -f "$file" ]; then
        echo "Error: The file '$file' does not exist."
        continue  # Skip to the next file
    fi

    # Count the number of lines in the file
    NumLines=$(wc -l < "$file")

    # Output the result
    echo "The file '$file' has $NumLines lines."

done

# Exit the script
exit 0
