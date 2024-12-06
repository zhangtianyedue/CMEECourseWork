"""
This script identifies oak trees from a list of species and demonstrates 
two methods for processing the data: using `for` loops and list comprehensions.

Features:
1. Identifies oak trees (genus `Quercus`) from a given list of taxa.
2. Demonstrates two approaches for:
   - Filtering oak species: 
     - Using `for` loops.
     - Using list comprehensions.
   - Converting oak species names to uppercase:
     - Using `for` loops.
     - Using list comprehensions.

Function:
- `is_an_oak(name)`: Checks whether a given species name belongs to the genus `Quercus`.

Purpose:
- To showcase basic data filtering and transformation using both `for` loops and list comprehensions.

Example Output:
- Oak species identified: {'Quercus robur', 'Quercus cerris', 'Quercus petraea'}
- Oak species in uppercase: {'QUERCUS ROBUR', 'QUERCUS CERRIS', 'QUERCUS PETRAEA'}
"""
## Finds just those taxa that are oak trees from a list of species

taxa = [ 'Quercus robur',
         'Fraxinus excelsior',
         'Pinus sylvestris',
         'Quercus cerris',
         'Quercus petraea',
       ]

def is_an_oak(name):
    return name.lower().startswith('quercus ')

##Using for loops
oaks_loops = set()
for species in taxa:
    if is_an_oak(species):
        oaks_loops.add(species)
print(oaks_loops)

##Using list comprehensions   
oaks_lc = set([species for species in taxa if is_an_oak(species)])
print(oaks_lc)

##Get names in UPPER CASE using for loops
oaks_loops = set()
for species in taxa:
    if is_an_oak(species):
        oaks_loops.add(species.upper())
print(oaks_loops)

##Get names in UPPER CASE using list comprehensions
oaks_lc = set([species.upper() for species in taxa if is_an_oak(species)])
print(oaks_lc)
