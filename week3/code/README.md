### 1. `oaks_debugme.py`
This script reads a CSV file of tree species and identifies which species belong to the Quercus genus (oak trees). It processes each row and filters out non-oak species, then writes the filtered data to an output file. 

- **Input**: `../data/TestOaksData.csv` (a CSV file containing genus and species information).
- **Output**: `../results/oaks_debugme_results.csv` (a CSV file listing all oak species from the input data).

The main functionality of this script is to determine if a given species belongs to the "Quercus" genus using the `is_an_oak` function, which checks whether the genus name starts with 'quercus'.

### 2. `align_seqs_fasta.py`
This script aligns two DNA sequences and finds the best possible alignment based on the number of matching bases. It reads the DNA sequences from FASTA files, aligns them at different starting points, and calculates a score for each alignment. The alignment with the highest score is then saved to an output file.

- **Input**: Two FASTA files containing DNA sequences, passed as command-line arguments or defaulted to `../data/407228326.fasta` and `../data/407228412.fasta`.
- **Output**: `../results/DNA_seq.txt` (a file containing the best alignment and its corresponding score).

The script uses the `calculate_score` function to align the sequences and compute the number of matches at various positions.

