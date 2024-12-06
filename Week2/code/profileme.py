"""
This script demonstrates the use of basic Python functions for list creation, string manipulation, and function execution.

Functions:
1. `my_squares(iters)`:
   - Creates a list of squares of numbers from 0 to `iters - 1` using a `for` loop.

2. `my_join(iters, string)`:
   - Concatenates a string `iters` times with a comma separator using a `for` loop.

3. `run_my_funcs(x, y)`:
   - Executes `my_squares` and `my_join` with the provided arguments and prints them.

Purpose:
- To illustrate iteration, string operations, and function composition in Python.
- Demonstrates potential inefficiencies in handling large-scale iterations and string concatenations.

Usage:
- Run the script to observe the operations, but note the potential for performance issues with large values of `x`.

"""

def my_squares(iters):  # Function to compute squares of numbers
    out = []  # Initialize an empty list
    for i in range(iters):  # Loop from 0 to iters - 1
        out.append(i ** 2)  # Append the square of the current number to the list
    return out  # Return the list of squares

def my_join(iters, string):  # Function to concatenate strings
    out = ''  # Initialize an empty string
    for i in range(iters):  # Loop from 0 to iters - 1
        out += string.join(", ")  # Join a comma separator and add it to the string
    return out  # Return the concatenated string

def run_my_funcs(x, y):  # Function to run both `my_squares` and `my_join`
    print(x, y)  # Print the input arguments for debugging
    my_squares(x)  # Call `my_squares` with argument `x`
    my_join(x, y)  # Call `my_join` with arguments `x` and `y`
    return 0  # Return 0 to indicate successful execution

run_my_funcs(10000000, "My string")  # Execute the functions with large input values
