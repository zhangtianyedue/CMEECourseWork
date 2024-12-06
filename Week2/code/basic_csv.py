"""
This script demonstrates reading from and writing to CSV files in Python using the `csv` module.

Features:
1. Reads a CSV file (`../data/testcsv.csv`) containing columns:
   - 'Species', 'Infraorder', 'Family', 'Distribution', 'Body mass male (Kg)'
2. Extracts and prints each row along with the species name.
3. Writes a new CSV file (`../data/bodymass.csv`) containing only the species names and their body masses.

Purpose:
- To illustrate how to handle CSV files, including reading data row by row and writing filtered content to a new file.
- Demonstrates the use of `csv.reader` for reading and `csv.writer` for writing.

Usage:
- Ensure the input file `../data/testcsv.csv` exists with the correct format.
- The script generates `../data/bodymass.csv` with two columns: species name and body mass.
"""
import csv

# Read a file containing:
# 'Species','Infraorder','Family','Distribution','Body mass male (Kg)'
with open('../data/testcsv.csv','r') as f:

    csvread = csv.reader(f)
    temp = []
    for row in csvread:
        temp.append(tuple(row))
        print(row)
        print("The species is", row[0])

# write a file containing only species name and Body mass
with open('../data/testcsv.csv','r') as f:
    with open('../data/bodymass.csv','w') as g:

        csvread = csv.reader(f)
        csvwrite = csv.writer(g)
        for row in csvread:
            print(row)
            csvwrite.writerow([row[0], row[4]])
