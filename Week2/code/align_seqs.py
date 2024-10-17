import csv

# Read the sequence from the CSV file
def read_sequences_from_csv(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        seq1 = ""#set up the empty string to recieve the new sequence 
        seq2 = ""
        for row in reader:
            if row['ID'] == 'seq1':
                seq1 = row['Sequence'].replace('=', '').replace('"', '').strip()  # "Remove the equal sign and the extra quotation marks." 
            elif row['ID'] == 'seq2':
                seq2 = row['Sequence'].replace('=', '').replace('"', '').strip()  # "Remove the equal sign and the extra quotation marks."
    return seq1, seq2


# caculate the score 
def calculate_score(s1, s2, l1, l2, startpoint):
    matched = ""  # Display the matching string.
    score = 0
    for i in range(l2):
        if (i + startpoint) < l1:
            if s1[i + startpoint] == s2[i]:  # If the bases match.
                matched = matched + "*"
                score = score + 1
            else:
                matched = matched + "-"
    return score, matched

# Find the best match and write the result to a CSV file.
def find_best_alignment(seq1, seq2, output_file):
    l1 = len(seq1)
    l2 = len(seq2)
    
    # if necessary, change the sequence 
    if l1 >= l2:
        s1 = seq1
        s2 = seq2
    else:
        s1 = seq2
        s2 = seq1
        l1, l2 = l2, l1  # swap the length of the sequence 
    
    my_best_align = None
    my_best_score = -1
    
    # find the perfect match 
    for i in range(l1):
        z, matched = calculate_score(s1, s2, l1, l2, i)
        if z > my_best_score:
            my_best_align = "." * i + s2  
            my_best_score = z

    # write the best result in the result file 
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Best Alignment", "Sequence 1", "Best Score"])
        writer.writerow([my_best_align, s1, my_best_score])

input_file = "../data/sequences.csv"
output_file = "../results/sequences_results.csv"

seq1, seq2 = read_sequences_from_csv(input_file)

print(f'seq1 = "{seq1}"')
print(f'seq2 = "{seq2}"')

# Find the best alignment and write the result to the output file.
find_best_alignment(seq1, seq2, output_file)

print("The result has been written:", output_file)


