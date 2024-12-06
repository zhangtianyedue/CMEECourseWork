#!/usr/bin/env python3

"""
lc1.py
This script demonstrates the use of list comprehensions and conventional loops
to extract information from a dataset of bird species. It creates three separate
lists: Latin names, common names, and mean body masses for each species.

Dataset:
- birds: A tuple of tuples containing species data (Latin name, common name, mean body mass).

Output:
- Prints the extracted lists for each approach.
"""

birds = ( 
    ('Passerculus sandwichensis', 'Savannah sparrow', 18.7),
    ('Delichon urbica', 'House martin', 19),
    ('Junco phaeonotus', 'Yellow-eyed junco', 19.5),
    ('Junco hyemalis', 'Dark-eyed junco', 19.6),
    ('Tachycineata bicolor', 'Tree swallow', 20.2),
)

def extract_list_comprehensions(data, index):
    """
    Extracts a specific field from the dataset using list comprehension.

    Parameters:
        data (tuple): Dataset containing bird species information.
        index (int): Index of the field to extract (0 for Latin name, 1 for common name, 2 for mean body mass).

    Returns:
        list: Extracted field values.
    """
    return [item[index] for item in data]

def extract_with_loops(data, index):
    """
    Extracts a specific field from the dataset using a conventional for loop.

    Parameters:
        data (tuple): Dataset containing bird species information.
        index (int): Index of the field to extract (0 for Latin name, 1 for common name, 2 for mean body mass).

    Returns:
        list: Extracted field values.
    """
    result = []
    for item in data:
        result.append(item[index])
    return result

def extract_with_while(data, index):
    """
    Extracts a specific field from the dataset using a while loop.

    Parameters:
        data (tuple): Dataset containing bird species information.
        index (int): Index of the field to extract (0 for Latin name, 1 for common name, 2 for mean body mass).

    Returns:
        list: Extracted field values.
    """
    result = []
    i = 0
    while i < len(data):
        result.append(data[i][index])
        i += 1
    return result

def main():
    """Main execution block."""
    # Using list comprehensions
    print("Step #1: Using List Comprehensions")
    print("Latin Names:", extract_list_comprehensions(birds, 0))
    print("Common Names:", extract_list_comprehensions(birds, 1))
    print("Mean Body Masses:", extract_list_comprehensions(birds, 2))
    print("-" * 40)

    # Using for loops
    print("Step #2: Using For Loops")
    print("Latin Names:", extract_with_loops(birds, 0))
    print("Common Names:", extract_with_loops(birds, 1))
    print("Mean Body Masses:", extract_with_loops(birds, 2))
    print("-" * 40)

    # Using while loops
    print("Step #3: Using While Loops")
    print("Latin Names:", extract_with_while(birds, 0))
    print("Common Names:", extract_with_while(birds, 1))
    print("Mean Body Masses:", extract_with_while(birds, 2))

if __name__ == "__main__":
    main()
