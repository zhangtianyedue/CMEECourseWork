"""
This script compares the performance of loops and list comprehensions for generating squares of numbers 
and the efficiency of different string concatenation methods. It uses functions from external modules 
`profileme` and `profileme2` for benchmarking.

Features:
1. Tests the performance of:
   - `my_squares` function implemented with loops (from `profileme`).
   - `my_squares` function implemented with list comprehensions (from `profileme2`).
2. Tests the performance of:
   - `my_join` function implemented with a loop-based join (from `profileme`).
   - `my_join` function implemented with the Python `join` method (from `profileme2`).

Purpose:
- To illustrate the differences in execution speed between loops and list comprehensions.
- To explore the efficiency of different string concatenation techniques.

Usage:
- Ensure the modules `profileme` and `profileme2` are available and include the corresponding functions.
- Use Python's `timeit` module to benchmark the performance of these functions with large inputs.
"""

iters = 1000000  # Number of iterations for the benchmarking tests

import timeit  # Import timeit for performance measurement

# Importing loop-based implementation of `my_squares` from `profileme`
from profileme import my_squares as my_squares_loops

# Importing list comprehension-based implementation of `my_squares` from `profileme2`
from profileme2 import my_squares as my_squares_lc

##############################################################################
# loops vs. the join method for strings: which is faster?
##############################################################################

mystring = "my string"  # String to be used for concatenation benchmarking

# Importing loop-based string join implementation from `profileme`
from profileme import my_join as my_join_join

# Importing string join implementation using Python's `join` method from `profileme2`
from profileme2 import my_join as my_join
