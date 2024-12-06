"""
This script demonstrates basic file input operations in Python, including 
reading a file line by line and handling blank lines.

Features:
1. Opens a file (`../sandbox/test.txt`) for reading.
2. Reads and prints all lines from the file using an implicit `for` loop.
3. Reopens the file to demonstrate skipping blank lines while printing only non-empty lines.

Purpose:
- To illustrate file handling in Python, including opening, reading, and closing files.
- Demonstrates how to process lines conditionally (e.g., ignoring blank lines).

Usage:
- Ensure the file `../sandbox/test.txt` exists with content to test the script.
"""
#############################
# FILE INPUT
#############################
# Open a file for reading
import os

# Define the path to the file
file_path = "../sandbox/test.txt"

# Check if the file exists
if not os.path.exists(file_path):
    # Ensure the sandbox directory exists
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    # Create the file and write content
    with open(file_path, 'w') as f:
        f.write("This is a test file.\n")
        f.write("Line 1: Hello, world!\n")
        f.write("Line 2: Python is awesome.\n")
    print(f"File created at: {file_path}")
else:
    print(f"File already exists at: {file_path}")

f = open('../sandbox/test.txt', 'r')
# use "implicit" for loop:
# if the object is a file, python will cycle over lines
for line in f:
    print(line)

# close the file
f.close()

# Same example, skip blank lines
f = open('../sandbox/test.txt', 'r')
for line in f:
    if len(line.strip()) > 0:
        print(line)

f.close()

"""
This script demonstrates basic file input operations in Python, including 
reading a file line by line and handling blank lines. ("with version")

Features:
1. Opens a file (`../sandbox/test.txt`) for reading.
2. Reads and prints all lines from the file using an implicit `for` loop.
3. Reopens the file to demonstrate skipping blank lines while printing only non-empty lines.

Purpose:
- To illustrate file handling in Python, including opening, reading, and closing files.
- Demonstrates how to process lines conditionally (e.g., ignoring blank lines).

Usage:
- Ensure the file `../sandbox/test.txt` exists with content to test the script.
""" 

#############################
# FILE INPUT
#############################
# Open a file for reading
with open('../sandbox/test.txt', 'r') as f:
    # use "implicit" for loop:
    # if the object is a file, python will cycle over lines
    for line in f:
        print(line)

# Once you drop out of the with, the file is automatically closed

# Same example, skip blank lines
with open('../sandbox/test.txt', 'r') as f:
    for line in f:
        if len(line.strip()) > 0:
            print(line)