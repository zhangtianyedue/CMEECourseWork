#!/usr/bin/env python3

"""A module with functions that demonstrate control flow, including sorting, comparing, and factorial calculations."""
__author__ = 'Zhang Tianye'

import sys

def foo_1(x):
    """
    Calculate the square root of a number.

    Parameters:
        x (float): A non-negative number.

    Returns:
        float: The square root of the number.
    """
    return x ** 0.5

def foo_2(x, y):
    """
    Return the larger of two numbers.

    Parameters:
        x (float): The first number.
        y (float): The second number.

    Returns:
        float: The larger of the two numbers.
    """
    return max(x, y)

def foo_3(x, y, z):
    """
    Return three numbers in ascending order.

    Parameters:
        x (float): The first number.
        y (float): The second number.
        z (float): The third number.

    Returns:
        list: A list of the three numbers in ascending order.
    """
    return sorted([x, y, z])

def factorial_recursive(n):
    """
    Calculate the factorial of a number using recursion.

    Parameters:
        n (int): A non-negative integer.

    Returns:
        int: The factorial of the number.
    """
    if n == 0 or n == 1:
        return 1
    return n * factorial_recursive(n - 1)

def factorial_iterative(n):
    """
    Calculate the factorial of a number using iteration.

    Parameters:
        n (int): A non-negative integer.

    Returns:
        int: The factorial of the number.
    """
    result = 1
    for i in range(2, n + 1):
        result *= i
    return result

def factorial_while_loop(n):
    """
    Calculate the factorial of a number using a while loop.

    Parameters:
        n (int): A non-negative integer.

    Returns:
        int: The factorial of the number.
    """
    result = 1
    while n > 1:
        result *= n
        n -= 1
    return result

def main(argv):
    """
    Test all the functions with example arguments.
    
    Parameters:
        argv (list): Command-line arguments (not used).

    Returns:
        int: Exit status (0 for success).
    """
    print(f"Square root of 25 is: {foo_1(25)}")
    print(f"Larger of 40 and 20 is: {foo_2(40, 20)}")
    print(f"Numbers 3, 2, 1 sorted are: {foo_3(3, 2, 1)}")
    print(f"Factorial of 5 using iteration is: {factorial_iterative(5)}")
    print(f"Factorial of 5 using recursion is: {factorial_recursive(5)}")
    print(f"Factorial of 5 using while loop is: {factorial_while_loop(5)}")
    return 0

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
