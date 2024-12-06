import csv
import os
print("Current working directory:", os.getcwd())

# Read the sequences from the CSV file
def read_sequences_from_csv(file_path):
    """
    Reads two sequences (seq1 and seq2) from a CSV file.

    Parameters:
        file_path (str): Path to the input CSV file.

    Returns:
        tuple: A tuple containing seq1 and seq2 as strings.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file does not contain the required sequences.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: Input file '{file_path}' not found.")

    seq1, seq2 = "", ""
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row['ID'] == 'seq1':
                seq1 = row['Sequence'].replace('=', '').replace('"', '').strip()
            elif row['ID'] == 'seq2':
                seq2 = row['Sequence'].replace('=', '').replace('"', '').strip()
    
    if not seq1 or not seq2:
        raise ValueError("Error: Input file must contain sequences with IDs 'seq1' and 'seq2'.")
    
    return seq1, seq2

# Calculate the alignment score
def calculate_score(s1, s2, l1, l2, startpoint):
    """
    Calculates the alignment score for two sequences.

    Parameters:
        s1 (str): The longer sequence.
        s2 (str): The shorter sequence.
        l1 (int): Length of s1.
        l2 (int): Length of s2.
        startpoint (int): Starting position for alignment.

    Returns:
        tuple: A tuple containing the score and the matched string.
    """
    matched = ""  # Display the matching string.
    score = 0
    for i in range(l2):
        if (i + startpoint) < l1:
            if s1[i + startpoint] == s2[i]:  # If the bases match.
                matched += "*"
                score += 1
            else:
                matched += "-"
    return score, matched

# Find the best alignment and write the result to a CSV file
def find_best_alignment(seq1, seq2, output_file):
    """
    Finds the best alignment between two sequences and writes the result to a CSV file.

    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        output_file (str): Path to the output CSV file.

    Returns:
        None
    """
    l1, l2 = len(seq1), len(seq2)
    
    # Ensure seq1 is the longer sequence
    if l1 < l2:
        seq1, seq2 = seq2, seq1
        l1, l2 = l2, l1  # Swap lengths
    
    best_align = None
    best_score = -1

    # Find the best alignment
    for i in range(l1):
        score, _ = calculate_score(seq1, seq2, l1, l2, i)
        if score > best_score:
            best_align = "." * i + seq2
            best_score = score

    # Write the result to the output file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Best Alignment", "Sequence 1", "Best Score"])
        writer.writerow([best_align, seq1, best_score])

# Main execution
if __name__ == "__main__":
    input_file = "../data/sequences.csv"
    output_file = "../results/sequences_results.csv"

    try:
        # Read sequences
        seq1, seq2 = read_sequences_from_csv(input_file)
        print(f'seq1 = "{seq1}"')
        print(f'seq2 = "{seq2}"')

        # Find the best alignment and write the result
        find_best_alignment(seq1, seq2, output_file)
        print("The result has been written:", output_file)
    except Exception as e:
        print(f"Error: {e}")



