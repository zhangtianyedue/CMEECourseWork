"""
Script to align two DNA sequences and find the best alignment based on matching bases.

This script reads two DNA sequences from FASTA files (or defaults to predefined file paths if none are provided),
compares them by aligning one sequence to the other at different starting points, and calculates an alignment score
based on the number of matching bases. The script identifies the best alignment (i.e., the alignment with the highest score)
and writes the result to a file.


Functions:
    - read_fasta: Reads a DNA sequence from a FASTA file, ignoring the header line.
    - calculate_score: Aligns two sequences starting from a specified position and returns the score based on matching bases.

Workflow:
    1. The script reads two DNA sequences from the command line or uses default files.
    2. It identifies the longer sequence and aligns the shorter one at various starting points.
    3. For each alignment, it computes a score based on how many bases match.
    4. It determines the alignment with the highest score and saves the result to an output file.

Input:
    - Two FASTA files containing DNA sequences, passed as command-line arguments or predefined file paths.
      (FASTA files should be in plain text format, without `.rtf` or other extensions.)

Output:
    - The best alignment and corresponding score are written to '../Results/DNA_seq.txt'.
    - Alignment details are printed to the terminal.

Author: Saskia Pearce (sp621@imperial.ac.uk)
Version: 3.9
"""

#input file name and seq1.fasta seq2.fasta 

import sys 
import os

# Define the base directory for your data (relative path)
base_dir = "../data/"

# Check if command-line arguments are provided
if len(sys.argv) >= 3: 
    # Use the paths provided by command-line arguments
    seq1  = os.path.join(base_dir, sys.argv[1])
    seq2  = os.path.join(base_dir, sys.argv[2])
else:
    # Use default paths if no arguments are provided
    seq1  = os.path.join(base_dir, "407228326.fasta")
    seq2  = os.path.join(base_dir, "407228412.fasta")

def read_fasta(filename): 
    with open(filename, 'r') as f: 
        lines = f.readlines()
    sequence = ''.join([line.strip() for line in lines if not line.startswith(">")])
    return sequence 

seq1 = read_fasta(seq1)
seq2 = read_fasta(seq2)

l1 = len(seq1) #assign longested sequence from file 
l2 = len(seq2) 
if l1 >= l2:
    s1 = seq1
    s2 = seq2
else:
    s1 = seq2
    s2 = seq1
    l1, l2 = l2, l1 # swap the two lengths

# A function that computes a score by returning the number of matches starting
# from arbitrary startpoint (chosen by user)
def calculate_score(s1, s2, l1, l2, startpoint):
    """
    Calculate the alignment score between two DNA sequences starting from a specified position.

    This function aligns a portion of sequence `s2` with sequence `s1` starting at the given 
    `startpoint`. It compares the bases of both sequences, adds to the score if they match, 
    and prints the alignment visually with '*' indicating a match and '-' indicating a mismatch.

    Args:
        s1 (str): The first DNA sequence (usually the longer one).
        s2 (str): The second DNA sequence (usually the shorter one to be aligned with `s1`).
        l1 (int): The length of sequence `s1`.
        l2 (int): The length of sequence `s2`.
        startpoint (int): The starting position in `s1` where the alignment with `s2` begins.

    Returns:
        int: The alignment score representing the number of matching bases between `s1` and `s2`.

    Example:
        >>> calculate_score("AGCTGAC", "GCT", 7, 3, 1)
        .***
        .GCT
        AGCTGAC
        3
    """
    matched = "" # to hold string displaying alignements
    score = 0
    for i in range(l2):
        if (i + startpoint) < l1:
            if s1[i + startpoint] == s2[i]: # if the bases match
                matched = matched + "*"
                score = score + 1
            else:
                matched = matched + "-"

    # some formatted output
    print("." * startpoint + matched)      #startpoitn has to shift up every time       
    print("." * startpoint + s2)
    print(s1)
    print(score) 
    print(" ")

    return score

# Test the function with some example starting points:
# calculate_score(s1, s2, l1, l2, 0)
# calculate_score(s1, s2, l1, l2, 1)
# calculate_score(s1, s2, l1, l2, 5)

my_best_align = None
my_best_score = -1

with open('../Results/DNA_seq.txt', 'w') as f:
    for i in range(l1):
        score = calculate_score(s1, s2, l1, l2, i)
        if score > my_best_score:
            my_best_align = "." * i + s2
            my_best_score = score
    f.write("Best alignment:\n")
    f.write(f"{my_best_align}\n{s1}\n")
    f.write(f"Best score: {my_best_score}\n")

print(my_best_align)
print(s1)
print("Best score:", my_best_score)




#results = '../results/alignment_results.txt'

# now try to find the best match (highest score) for the two sequences
#my_best_align = None
#my_best_score = -1 #alignment will always be higher than this starting score

#f = open('../Results/DNA_seq.txt','w')
#for i in range(l1): # Note that you just take the last alignment with the highest score
    #z = calculate_score(s1, s2, l1, l2, i) #arguement in the file 
    #if z > my_best_score:
        #my_best_align = "." * i + s2 # think about what this is doing!
        #my_best_score = z 
        #f.write("Best alignment:\n")
        #f.write(f"{my_best_align}\n{s1}\n")
        #f.write(f"Best score: {my_best_score}\n")

#f.close()

#print(my_best_align)
#print(s1)
#print("Best score:", my_best_score) # best score printed

# assigning seq1.fasta (first argument)






# assigning seq2.fasta (second arguement)



