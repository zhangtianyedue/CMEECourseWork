import csv
import sys
from fuzzywuzzy import fuzz
#because im using the Macos system,to bypass the limitations of Homebrew, I create a virtual environment to install the required dependencies. The steps are as follows:
#1.Create a virtual environment:
#bash
#python3 -m venv myenv
##2.Activate the virtual environment:
#bash
#source myenv/bin/activate
##3.Install the required module:
#bash
#pip install fuzzywuzzy
#4.Run my script:
#python3 oaks_debugme.py
#5.exit the virtual enviroment 
#deactivate


"""
oaks_debugme.py
This script processes a CSV file of tree species and identifies entries that are oak trees based on fuzzy string matching.
It writes the oak tree entries to a new CSV file.

Dependencies:
- fuzzywuzzy (install via pip: pip install fuzzywuzzy)

Usage:
Ensure the input file `TestOaksData.csv` is in the ../data/ directory, and the output will be saved to ../results/oak_result.csv.
"""

# Define function
def is_an_oak(name):
    """
    Returns True if name starts with 'quercus'
    
    >>> is_an_oak('Quercus robur')
    True
    >>> is_an_oak('Fagus sylvatica')
    False
    >>> is_an_oak('Quercuss')
    True
    """
    return fuzz.partial_ratio(name.lower(), 'quercus') > 85  

def main(argv):
    """
    Main function to process the tree species data and write oak entries to a new file.
    """
    try:
        # Open input and output files
        f = open('../data/TestOaksData.csv', 'r')
        g = open('../results/oak_result.csv', 'w', newline='')
    except FileNotFoundError as e:
        print(f"Error: {e}. Please check if the file paths are correct.")
        return 1
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return 1
    
    taxa = csv.reader(f)
    csvwrite = csv.writer(g)
    oaks = set()
    try:
        for row in taxa:
            print(row)
            print("The genus is: ")
            print(row[0] + '\n')
            if is_an_oak(row[0]):
                print('FOUND AN OAK!\n')
                csvwrite.writerow([row[0], row[1]])
    except IndexError:
        print("Error: The input file may have an incorrect format. Ensure it contains at least two columns.")
        return 1
    except Exception as e:
        print(f"An error occurred while processing the file: {e}")
        return 1
    finally:
        # Close the files
        f.close()
        g.close()

    return 0

if __name__ == "__main__":
    status = main(sys.argv)
