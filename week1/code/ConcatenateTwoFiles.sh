#!/bin/bash

# Original functionality: merging two files into a third one

# Check if three arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 file1 file2 output_file"
    exit 1
fi

# Check if the first input file exists
if [ ! -f "$1" ]; then
    echo "Error: File $1 does not exist."
    exit 1
fi

# Check if the second input file exists
if [ ! -f "$2" ]; then
    echo "Error: File $2 does not exist."
    exit 1
fi

# Check if the output file exists
if [ -f "$3" ]; then
    # If output file exists, ask user if they want to overwrite it
    echo "Warning: Output file $3 already exists. Do you want to overwrite it? (y/n)"
    read answer
    if [ "$answer" != "y" ]; then
        echo "Exiting without overwriting the file."
        exit 0
    fi
else
    # If the output file does not exist, create it
    touch "$3"
    echo "Output file $3 does not exist. Creating it..."
fi

# Merge the first input file into the output file (overwrite)
cat "$1" > "$3"

# Append the second input file to the output file
cat "$2" >> "$3"

# Display the result
echo "Merged File is:"
cat "$3"

exit 0


