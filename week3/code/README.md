***All the codes this week use relative paths. Therefore, please make sure to set the working directory to `week3/code` when running the scripts. Otherwise, the scripts may fail to execute due to missing files.***

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

***NEW VERSION update on 03/12/2024
The above explains the group assignment. In addition, Week 3 also includes class learning materials related to R. Below are some updates.

For the `browse.R` file, special attention is needed because the `browser()` function is used. When running the code, the program will pause when it reaches `browser()`. 
In the terminal, it will display `Browse[2]>`, indicating that you are debugging the `Exponential` function. In debug mode, you can inspect the variable values within the function. For example, typing `N` will allow you to check the current state of `N`.

The script, `DataWrang.R`, `DataWrangTidy.R`requires the `tidyverse` package for data manipulation and analysis. Please ensure that the `tidyverse` package is installed in your R environment before running the script. If it is not already installed, you can install it by running the following command in your R console:
```r
install.packages("tidyverse")
``` Once installed, the package will be loaded automatically when the script is executed.