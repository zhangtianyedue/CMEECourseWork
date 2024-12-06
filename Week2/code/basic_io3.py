"""
This script demonstrates how to store and retrieve Python objects 
using the `pickle` module for serialization and deserialization.

Features:
1. Creates a dictionary object with key-value pairs.
2. Saves the dictionary to a binary file (`../sandbox/testp.p`) using `pickle.dump`.
3. Loads the dictionary back into a new variable using `pickle.load`.
4. Prints the retrieved dictionary to verify the stored data.

Purpose:
- To illustrate object serialization (saving Python objects to files) and deserialization (loading objects back into memory).
- Demonstrates how to handle binary file operations with the `pickle` module.

Usage:
- The script saves the dictionary object to `../sandbox/testp.p` and reloads it for verification.
- Ensure the `../sandbox/` directory exists to avoid file errors.
"""
#############################
# STORING OBJECTS
#############################
# To save an object (even complex) for later use
my_dictionary = {"a key": 10, "another key": 11}

import pickle

f = open('../sandbox/testp.p','wb') ## note the b: accept binary files
pickle.dump(my_dictionary, f)
f.close()

## Load the data again
f = open('../sandbox/testp.p','rb')
another_dictionary = pickle.load(f)
f.close()

print(another_dictionary)
