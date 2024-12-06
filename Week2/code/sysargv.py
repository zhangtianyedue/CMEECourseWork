"""
This script demonstrates how to access command-line arguments in Python using the `sys` module.

It prints:
1. The name of the script being executed.
2. The number of arguments passed to the script.
3. The list of all arguments as a string.

Usage:
- Run this script with additional arguments to see how they are handled.
"""
#!/usr/bin/env python3

import sys
print("This is the name of the script: ", sys.argv[0])
print("Number of arguments: ", len(sys.argv))
print("The arguments are: " , str(sys.argv))