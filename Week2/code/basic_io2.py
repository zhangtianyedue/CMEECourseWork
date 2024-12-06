"""
This script demonstrates basic file output operations in Python, including 
saving the elements of a list to a file.

Features:
1. Creates a list of numbers from 0 to 99 using `range(100)`.
2. Opens a file (`../sandbox/testout.txt`) for writing.
3. Writes each element of the list to the file, adding a newline character after each number.
4. Closes the file to ensure the data is saved properly.

Purpose:
- To illustrate how to write data to a file in Python.
- Demonstrates handling line-by-line file writing and proper file closure.

Usage:
- The script generates a file `../sandbox/testout.txt` containing numbers from 0 to 99, each on a separate line.
"""
#############################
# FILE OUTPUT
#############################
# Save the elements of a list to a file
list_to_save = range(100)

f = open('../sandbox/testout.txt','w')
for i in list_to_save:
    f.write(str(i) + '\n') ## Add a new line at the end

f.close()