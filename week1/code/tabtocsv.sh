#!/bin/sh
# Author: Your name you.login@imperial.ac.uk
# Script: tabtocsv.sh
# Description: substitute the tabs in the files with commas
# Saves the output into a .csv file
# Arguments: 1 -> tab delimited file
# Date: Oct 2024

# Check if an argument (file) is provided
if [ $# -eq 0 ]; then
    # If no file is provided, print an error message and exit
    echo "Error: No file provided."
    echo "Usage: ./tabtocsv.sh <file>"
    exit 1  # Exit with a non-zero status to indicate an error
fi

# Check if the input file has a .csv extension
if [ "${1##*.}" != "csv" ]; then
    # If the file does not end with .csv, print an error and exit
    echo "Error: The file is not a .csv file."
    exit 1  # Exit with a non-zero status to indicate an error
fi

# If checks are passed, execute the original logic
echo "Creating a comma delimited version of $1 ..."
# Use 'cat' to read the file, 'tr' to replace tabs with commas, and append the result to a new file
cat $1 | tr -s "\t" "," >> "${1%.csv}_converted.csv"
echo "Done!"
exit 0  # Exit successfully

