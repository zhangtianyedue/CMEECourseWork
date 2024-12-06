#!/usr/bin/env python3

"""
This script processes a tuple of bird species, where each entry contains the Latin name,
common name, and body mass of a bird. It prints each species' information on a separate line
in a formatted manner.

Dataset:
- birds: A tuple of tuples containing (Latin name, common name, mass).

Output:
- Prints each bird's details: Latin name, common name, and mass.
"""

birds = (
    ('Passerculus sandwichensis', 'Savannah sparrow', 18.7),
    ('Delichon urbica', 'House martin', 19),
    ('Junco phaeonotus', 'Yellow-eyed junco', 19.5),
    ('Junco hyemalis', 'Dark-eyed junco', 19.6),
    ('Tachycineata bicolor', 'Tree swallow', 20.2),
)

for i in birds:
    print("Latin name:", i[0], "Common name:", i[1], "Mass:", i[2])
