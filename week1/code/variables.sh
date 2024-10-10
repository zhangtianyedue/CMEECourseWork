#!/bin/sh

# This script demonstrates the use of special variables, assigned variables, reading input, and handling conditions.

# Check if there are any parameters passed to the script
if [ $# -eq 0 ]; then
    echo "No parameters were passed to the script."
else
    # Display the number of parameters passed to the script
    echo "This script was called with $# parameters."
fi

# Display the name of the script
echo "The script's name is $0"

# If there are any arguments, display them
if [ $# -gt 0 ]; then
    echo "The arguments are: $@"
fi

# Check if at least two arguments were passed
if [ $# -ge 2 ]; then
    # Display the first and second arguments
    echo "The first argument is $1"
    echo "The second argument is $2"
else
    echo "Not enough arguments were passed to display the first and second arguments."
fi

# Assigned Variables: Explicit declaration of a string variable
MY_VAR='some string'
echo 'The current value of the variable is:' $MY_VAR
echo

# Prompt the user to input a new value for the variable
echo 'Please enter a new string:'
read MY_VAR

# Display the new value entered by the user
echo
echo 'The new value of the variable is:' $MY_VAR
echo

# Handling multiple values from user input
echo 'Please enter two numbers separated by space(s):'
read a b

# Check if both inputs are provided
if [ -z "$a" ] || [ -z "$b" ]; then
    echo "You did not provide two numbers."
else
    # Check if the inputs are valid numbers
    if ! [[ "$a" =~ ^[0-9]+$ ]] || ! [[ "$b" =~ ^[0-9]+$ ]]; then
        echo "Error: Both inputs must be numbers."
    else
        # Display the two numbers entered by the user
        echo 'You entered' $a 'and' $b '.'

        # Perform arithmetic using command substitution
        MY_SUM=$(expr $a + $b)
        echo "Their sum is: $MY_SUM"
    fi
fi

exit 0
