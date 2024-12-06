"""
A simple script to demonstrate the behavior of the `__name__` variable in Python.

When this script is executed directly, it identifies itself as the main program 
and prints a specific message. If it is imported as a module into another program, 
it acknowledges that it is being used as an import. Additionally, it prints the 
current value of the `__name__` variable in both cases.

Usage:
- Run the script directly to see how `__name__` behaves.
- Import it into another script to observe the difference in output.
"""
#!/usr/bin/env python3
# Filename: using_name.py

if __name__ == '__main__':
    print('This program is being run by itself!')
else:
    print('I am being imported from another script/program/module!')

print("This module's name is: " + __name__)