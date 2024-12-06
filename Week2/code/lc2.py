#!/usr/bin/env python3

"""
lc2.py
This script demonstrates the use of list comprehensions and conventional loops
to filter UK rainfall data by specific conditions (e.g., rainfall > 100mm or < 50mm).

Dataset:
- rainfall: A tuple of tuples containing monthly rainfall data (month, rainfall in mm).

Output:
- Lists of months and rainfall values meeting specified conditions.
"""

rainfall = (
    ('JAN', 111.4), ('FEB', 126.1), ('MAR', 49.9), ('APR', 95.3),
    ('MAY', 71.8), ('JUN', 70.2), ('JUL', 97.1), ('AUG', 140.2),
    ('SEP', 27.0), ('OCT', 89.4), ('NOV', 128.4), ('DEC', 142.2),
)

def filter_rainfall(data, condition):
    """
    Filters rainfall data based on a condition.

    Parameters:
        data (tuple): A tuple of tuples containing rainfall data (month, rainfall).
        condition (function): A function defining the filter condition.

    Returns:
        list: A list of tuples or month names meeting the condition.
    """
    return [item for item in data if condition(item)]

def filter_rainfall_loops(data, condition):
    """
    Filters rainfall data using a for loop based on a condition.

    Parameters:
        data (tuple): A tuple of tuples containing rainfall data (month, rainfall).
        condition (function): A function defining the filter condition.

    Returns:
        list: A list of tuples or month names meeting the condition.
    """
    result = []
    for item in data:
        if condition(item):
            result.append(item)
    return result

def main():
    """Main execution block."""
    # List comprehension examples
    print("Step #1: Using List Comprehensions")
    print("Months and rainfall > 100mm:", filter_rainfall(rainfall, lambda x: x[1] > 100))
    print("Months with rainfall < 50mm:", [month[0] for month in filter_rainfall(rainfall, lambda x: x[1] < 50)])
    print("-" * 40)

    # Conventional loops
    print("Step #2: Using For Loops")
    print("Months and rainfall > 100mm:", filter_rainfall_loops(rainfall, lambda x: x[1] > 100))
    print("Months with rainfall < 50mm:", [month[0] for month in filter_rainfall_loops(rainfall, lambda x: x[1] < 50)])
    print("-" * 40)

if __name__ == "__main__":
    main()
