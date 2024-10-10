#!/bin/bash

# Ensure one argument is provided (the input file)
if [ $# -ne 1 ]; then
    echo "False, Please input the inputfilename"
    exit 1
fi

input_file="$1"
output_dir="/Users/tianyezhang/Desktop/CMEECourseWork/week1/results"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "False, the inputfile $input_file does not exist."
    exit 1
fi

# Extract the base name of the input file (without the extension)
base_name=$(basename "$input_file" .csv)

# Automatically set the output file name in the results directory
output_file="$output_dir/${base_name}result.csv"

# Ensure the results directory exists; create it if it doesn't
mkdir -p "$output_dir"

# Check if the file is already space-separated (no commas)
if ! grep -q ',' "$input_file"; then
    echo "The file is already space-separated. Copying the file to the output."
    cp "$input_file" "$output_file"
    echo "File copied successfully: $output_file"
else
    # Replace commas with spaces and save to output file
    tr ',' ' ' < "$input_file" > "$output_file"
    echo "Transformation successful and saved as: $output_file"
fi

exit 0










