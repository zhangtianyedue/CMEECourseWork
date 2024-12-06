"""
This script defines and demonstrates the use of a simple function `foo`.

Functionality:
- The function `foo` takes a single argument `x`, squares it (`x = x * x`), and prints the result.

Purpose:
- To demonstrate basic function definition, parameter manipulation, and printing in Python.

Example:
- Calling `foo(2)` prints `4` as the result.
"""

def foo(x):
    x *= x # same as x = x*x
    print(x)

foo(2)