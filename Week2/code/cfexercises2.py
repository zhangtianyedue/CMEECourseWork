"""
This script defines and demonstrates various `hello_*` functions, each showcasing different 
control flow techniques in Python, including `for` loops, `while` loops, conditional statements, 
and the use of `break`.

Functions:
1. `hello_1(x)`: Prints 'hello' for every number from 0 to `x-1` divisible by 3.
2. `hello_2(x)`: Prints 'hello' for every number from 0 to `x-1` where:
   - The remainder when divided by 5 is 3, or
   - The remainder when divided by 4 is 3.
3. `hello_3(x, y)`: Prints 'hello' for every number from `x` to `y-1`.
4. `hello_4(x)`: Prints 'hello' repeatedly, incrementing `x` by 3 until `x` equals 15.
5. `hello_5(x)`: Prints 'hello' under specific conditions while incrementing `x`:
   - Prints 'hello' 7 times if `x` equals 31.
   - Prints 'hello' once if `x` equals 18.
6. `hello_6(x, y)`: Demonstrates a `while` loop with a `break` condition:
   - Continuously prints 'hello' with the value of `y`, incrementing `y` by 1, 
     until `y` reaches 6.

Purpose:
- To illustrate various loop constructs and conditional branching in Python.
- Provides examples of `for` and `while` loops, nested loops, and breaking out of loops.

Example:
- Each function is called with sample inputs, and the corresponding output is printed to the console.
"""

########################
def hello_1(x):
    for j in range(x):
        if j % 3 == 0:
            print('hello')
    print(' ')

hello_1(12)

########################
def hello_2(x):
    for j in range(x):
        if j % 5 == 3:
            print('hello')
        elif j % 4 == 3:
            print('hello')
    print(' ')

hello_2(12)

########################
def hello_3(x, y):
    for i in range(x, y):
        print('hello')
    print(' ')

hello_3(3, 17)

########################
def hello_4(x):
    while x != 15:
        print('hello')
        x = x + 3
    print(' ')

hello_4(0)

########################
def hello_5(x):
    while x < 100:
        if x == 31:
            for k in range(7):
                print('hello')
        elif x == 18:
            print('hello')
        x = x + 1
    print(' ')

hello_5(12)

# WHILE loop with BREAK
def hello_6(x, y):
    while x: # while x is True
        print("hello! " + str(y))
        y += 1 # increment y by 1 
        if y == 6:
            break
    print(' ')

hello_6 (True, 0)