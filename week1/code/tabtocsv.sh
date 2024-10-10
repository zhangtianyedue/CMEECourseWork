#!/bin/sh
# Author: zhangtianye
# Script: tabtocsv.sh
# Description: Substitute the tabs in the files with commas
#
# Saves the output into a .csv file
# Arguments: 1 -> tab delimited file
# Date: Oct 2024

# Check if the file already exist
if [ $# -lt 1 ]; then
    echo "Error: Please enter an input file name."
    exit 1
fi
echo "Creating a comma delimited version of $1 ..."
cat "$1" | tr -s "\t" "," >> "$1.csv"
echo "Done!"
exit
