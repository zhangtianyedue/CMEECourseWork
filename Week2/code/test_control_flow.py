#!/usr/bin/env python3

"""
This script demonstrates the use of control statements in Python through 
a function that determines whether a number is even or odd. It includes 
embedded tests using the `doctest` module.

Features:
1. `even_or_odd(x)`: A function to determine if a number is even or odd, 
   with examples tested using `doctest`.
2. `main(argv)`: A function to execute and print results for example inputs.
3. Integration with `doctest` to validate function behavior against embedded examples.

Usage:
- Run the script to see results for sample inputs.
- The `doctest` module will automatically test the embedded examples.
"""


__author__ = 'zhangtianye'
__version__ = '0.0.1'

import sys
import doctest # Import the doctest module

def even_or_odd(x=0):
    """Find whether a number x is even or odd.
      
    >>> even_or_odd(10)
    '10 is Even!'
    
    >>> even_or_odd(5)
    '5 is Odd!'
        
    in case of negative numbers, the positive is taken:    
    >>> even_or_odd(-2)
    '-2 is Even!'
    
    """
    #Define function to be tested
    if x % 2 == 0:
        return f"{x} is Even!"
    return f"{x} is Odd!"

def main(argv): 
    print(even_or_odd(22))
    print(even_or_odd(33))
    return 0

if (__name__ == "__main__"):
    status = main(sys.argv)

doctest.testmod()   # To run with embedded tests