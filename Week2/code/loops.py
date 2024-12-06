# FOR loops
"""
This script demonstrates the use of basic control flow structures in Python: `for` and `while` loops.

Features:
1. `FOR` loops:
   - Iterates over a range of numbers (0 to 4) and prints each value.
   - Iterates over a list containing mixed data types and prints each element.
   - Iterates over a list of numbers, calculates their cumulative sum, and prints the total at each step.

2. `WHILE` loop:
   - Iteratively increments a variable until it reaches 100, printing each value.

Purpose:
- To illustrate how loops operate in Python with practical examples.
"""

for i in range(5):
    print(i)
my_list = [0, 2, "geronimo!", 3.0, True, False]
for k in my_list:
    print(k)

total = 0
summands = [0, 1, 11, 111, 1111]
for s in summands:
    total = total + s
    print(total)

# WHILE loop
z = 0
while z < 100:
    z = z + 1
    print(z)

    