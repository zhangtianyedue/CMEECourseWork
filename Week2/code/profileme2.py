"""
This script demonstrates optimized Python functions for generating lists of squares and concatenating strings.
It uses a list comprehension for efficiency in square generation and a loop for string concatenation.

Functions:
1. `my_squares(iters)`:
   - Creates a list of squares of numbers from 0 to `iters - 1` using a list comprehension.

2. `my_join(iters, string)`:
   - Concatenates the given string `iters` times, separated by commas, using a `for` loop.

3. `run_my_funcs(x, y)`:
   - Calls `my_squares` and `my_join` with the provided arguments and prints the inputs for debugging.

Purpose:
- To illustrate basic iteration techniques for list generation and string manipulation.
- Highlights the trade-offs between efficiency and clarity in Python loops.

Usage:
- Run the script to observe the behavior of the functions with large inputs.

Note:
- Be cautious of performance limitations for large values of `iters` due to memory usage and processing time.
"""

def my_squares(iters):  # Function to compute squares of numbers
    out = [i ** 2 for i in range(iters)]  # Use a list comprehension to generate squares from 0 to iters - 1
    return out  # Return the list of squares

def my_join(iters, string):  # Function to concatenate a string multiple times
    out = ''  # Initialize an empty string
    for i in range(iters):  # Loop through the range of iterations
        out += ", " + string  # Append a comma and the given string to the output
    return out  # Return the concatenated string

def run_my_funcs(x, y):  # Function to execute both `my_squares` and `my_join`
    print(x, y)  # Print the input arguments for debugging purposes
    my_squares(x)  # Call `my_squares` with argument `x`
    my_join(x, y)  # Call `my_join` with arguments `x` and `y`
    return 0  # Return 0 to indicate successful execution

run_my_funcs(10000000, "My string")  # Execute the functions with a large number of iterations and a test string
