#The 'code' folder contains a total of 7 Python files.


***For this assignment, I have strictly used relative paths. Therefore, please ensure that when running any files, your working directory is set to `CMEECourseWork/Week2/code/`. If this is not done, the scripts will fail to locate the required files and display a "file not found" error.***


update version on 01/12/2024
Among them, the code for the 1st question **lc1.py** demonstrates the extraction of Latin names, common names, and mean body masses of bird species using list comprehensions, `for` loops, and `while` loops. The code has been optimized by modularizing repetitive logic into parameterized functions (`extract_list_comprehensions`, `extract_with_loops`, and `extract_with_while`), reducing redundancy and improving maintainability. Clear outputs and structured sections highlight each method, making the script both educational and efficient.

The code for the 1st question **lc1.py** filters UK rainfall data from 1910 to find months with rainfall > 100mm and < 50mm. It demonstrates two methods: `for` loops and `while` loops, showcasing both approaches to complete the task efficiently and clearly. The results are printed for easy comparison.

The code for the 1st question **tuple.py** processes a dataset of bird species, where each entry includes the Latin name, common name, and body mass. It formats and prints each bird's details (Latin name, common name, and mass) on a separate line for clarity. The data is stored in a tuple of tuples and iterated over using a simple `for` loop.

The code for the second question is **cfexercises1.py**, which demonstrates control flow through various functions, including square root calculation, finding the maximum of two numbers, sorting three numbers, and calculating factorials using recursive, iterative, and while-loop methods. The improvements include reducing redundancy by refactoring factorial functions (`foo_4`, `foo_5`, `foo_6`) into three distinct, well-named methods (`factorial_iterative`, `factorial_recursive`, `factorial_while_loop`) while retaining their unique logic. Built-in functions like `max` and `sorted` were used to simplify `foo_2` and `foo_3`, improving readability and efficiency. Detailed docstrings were also added for clarity and maintainability.

The code for the first question **dictionary.py** creates a dictionary that maps order names (categories) to sets of species using a provided list of tuples as input. Each tuple contains a species name and its corresponding category. The script processes this data using two methods: one with defaultdict for automatic set initialization and another with a dictionary comprehension for a concise alternative. The input data is hardcoded in the taxa list, and the output is a dictionary printed directly to the console, showing the mapping of categories to species.


The code for the third question is named** align_seqs.py**. Please note that when running this file, make sure you are in the **week2/code directory**, as I used relative paths when importing files. The file being imported is named **sequences.csv**, which is located in the**week2/data** folder. The final result will be directly written to the week2/results folder, with the filename **sequences_results.csv**，This script performs sequence alignment for two DNA sequences provided in a CSV file (../data/sequences.csv) and writes the best alignment and score to ../results/sequences_results.csv. It calculates the alignment using a sliding window approach and ensures that missing input files or formatting issues are handled gracefully. The input file must have seq1 and seq2 in the ID column, with their sequences in the Sequence column. Run the script from the code directory using python sequence_alignment.py. The results folder will be created automatically if it does not exist.


The code for the fourth question is named oaks_debugme.py. Please note that when running this file, make sure you are in the week2/code directory, as I used relative paths when importing files. The file being imported is named TestOaksData.csv, which is located in the week2/data folder. The final result will be directly written to the week2/results folder, with the filename oak_result.csv.
Dependencies: This script requires the fuzzywuzzy package for fuzzy string matching. You can install it using the following command:
bash
pip install fuzzywuzzy

because im using the Macos system,to bypass the limitations of Homebrew, I create a virtual environment to install the required dependencies. The steps are as follows:
1.Create a virtual environment:
bash
python3 -m venv myenv
2.Activate the virtual environment:
bash
source myenv/bin/activate
3.Install the required module:
bash
pip install fuzzywuzzy
4.Run my script:
python3 oaks_debugme.py
5.exit the virtual enviroment 
deactivate




update version on  06/12/2024
This folder contains all the Python files that were introduced as examples during class exercises. I have included them in this collection and taken the time to carefully understand each script. For all files, I have added script-level docstrings to clearly explain their purpose, functionality, and usage. These docstrings aim to make the code more accessible and understandable for anyone reviewing or using the scripts. This effort is part of my commitment to improving the readability, clarity, and educational value of these examples.

