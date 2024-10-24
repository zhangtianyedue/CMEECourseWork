"""
This module processes a CSV file of tree species, checks whether each species belongs to the 'Quercus' genus (oak trees),
and writes the oak species to an output CSV file, while filtering out duplicate entries.

Functions:
    - is_an_oak(name): Determines if a species belongs to the Quercus genus (oak trees).
    - main(argv): Reads input CSV data, processes the species list, and writes oaks to the results file.
    
Module-Level:
    The script performs the following tasks:
    1. Reads a CSV file ('../data/TestOaksData.csv') containing tree species data with 'Genus' and 'Species' columns.
    2. Identifies oak species (those from the 'Quercus' genus) and writes them into a new CSV file ('../results/oaks_debugme_results.csv').
    3. Ensures there are no duplicate oak species in the output.
    4. Handles both the presence and absence of headers in the input file.

    The script also includes unit tests within the `is_an_oak` function, which can be executed via Python's doctest module.

Usage:
    Run the script from the command line or from within another Python program.
    Example:
        $ python oaks_debugme.py
"""

import csv
import sys
import doctest

#==================================================================================================

# Define the function `is_an_oak`
def is_an_oak(name):
    """
    Determines if a tree species belongs to the Quercus genus (oak trees).
    
    Args:
        name (str): The name of the genus to check (usually the genus of a tree).
    
    Returns:
        bool: True if the species belongs to the Quercus genus, otherwise False.
    
    This function converts the input string to lowercase, checks if the name starts with 'quercus', and ensures
    that either the name exactly matches 'quercus' or the next character after 'quercus' is a space.
    
    Examples:
        >>> is_an_oak('Fagus sylvatica')
        False
        >>> is_an_oak('Quercus robur')
        True
        >>> is_an_oak('quercus ROBUR')
        True
        >>> is_an_oak('Quercuss robur')
        False
        >>> is_an_oak('Quercus')
        True
        >>> is_an_oak(' Quercus robur')
        False
    """
    # Convert the name to lowercase for case-insensitive matching
    name = name.lower()
    # Check if the name starts with 'quercus' and either matches exactly or has a space after 'quercus'
    return name.startswith('quercus') and (len(name) == len('quercus') or name[len('quercus')] == ' ')

#==================================================================================================

def main(argv):
    """
    Main function that processes a CSV file containing genus and species data, identifies oak species,
    and writes the results to an output CSV file.
    
    Args:
        argv (list): Command line arguments passed to the script (though not used directly in this implementation).
    
    Steps:
        1. Reads the input CSV file ('../data/TestOaksData.csv') containing tree species information.
        2. Identifies species that belong to the Quercus genus using the `is_an_oak` function.
        3. Writes oak species into the output file ('../results/oaks_debugme_results.csv'), avoiding duplicates.
        4. Handles input files with or without headers. If headers exist, they are written to the output.

    Input:
        The input CSV file should have two columns: 'Genus' and 'Species'.
    
    Output:
        The results CSV file will contain oak species filtered by genus and species.

    Returns:
        int: Status code (0 if successful).
    
    Example usage:
        $ python oaks_debugme.py
    
    Notes:
        - The input file '../data/TestOaksData.csv' is expected to exist and be properly formatted.
        - This function handles duplicates and only outputs unique oak species.
    """
    # Open the input CSV file with tree species data and prepare to write results to another CSV file
    f = open('../data/TestOaksData.csv', 'r')  # Input file
    g = open('../results/oaks_debugme_results.csv', 'w')  # Output file
    taxa = csv.reader(f)  # Read the input file
    csvwrite = csv.writer(g)  # Prepare to write to the output file
    oaks = set()  # Set to store unique oak species names

    # Skip the first row (header), if present, and write it to the output file
    header = next(taxa)
    if header[0].strip().lower() == 'genus' and header[1].strip().lower() == 'species':
        print("Skipping header row...")  # Inform that the header row is being skipped
        csvwrite.writerow(header)  # Write the header to the output file
    else:
        taxa = [header] + list(taxa)  # If no header row, add it back to the remaining rows
        header = ["Genus", "Species"]
        csvwrite.writerow(header)  # Write default header row

    # Iterate over each row and check if it belongs to the Quercus genus (oak trees)
    for row in taxa:
        print(row)
        print("The genus is: ", row[0])

        # If it's an oak, and not a duplicate, write it to the results file
        if is_an_oak(row[0]):
            rowfullname = row[0] + " " + row[1]
            if rowfullname not in oaks:  # Prevent duplicates
                print(rowfullname, " is an oak! \n")
                csvwrite.writerow([row[0], row[1]])  # Write the row to the results file
                oaks.add(rowfullname)  # Add the oak species to the set
        else:
            rowfullname = row[0] + " " + row[1]
            print(rowfullname, "is not an oak! \n")  # If not an oak, display message and skip

    return 0  # Return 0 to indicate successful execution

#==================================================================================================

if __name__ == "__main__":
    """
    Calls the main function and runs the doctest module to test the is_an_oak function.
    The script will exit with status code 0 upon successful execution.
    """
    import doctest
    doctest.testmod()  # Run doctests to validate the is_an_oak function
    status = main(sys.argv)  # Run the main function
    sys.exit(status)  # Exit the program with the returned status





