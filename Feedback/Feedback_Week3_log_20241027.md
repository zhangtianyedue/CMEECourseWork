
# Feedback on Project Structure, Workflow, and Code Structure

**Student:** Tianye Zhang

---

## General Project Structure and Workflow

- **Directory Organization**: The project structure is well-organized with clear weekly folders (`week1`, `Week2`, `week3`) and subdirectories (`code`, `data`, `results`, `sandbox`) under `week3`. This setup promotes easy navigation and separates data, code, and outputs effectively.
- **README Files**: The main README file (`README.log`) is present, although it lacks content. In `week3`, the `README.md` file is informative, detailing script functionality, input/output paths, and usage. Adding more details such as usage examples with sample inputs/outputs would improve accessibility.

### Suggested Improvements:
1. **Expand README Files**: Provide examples of input/output and usage for scripts, particularly `oaks_debugme.py` and `align_seqs_fasta.py`.
2. **Add `.gitignore`**: Including a `.gitignore` file to exclude temporary files and `.DS_Store` files would keep the repository clean.

## Code Structure and Syntax Feedback

### Python Scripts in `week3/code`

1. **align_seqs_fasta.py**:
   - **Purpose**: This script aligns two DNA sequences, determining the best alignment based on matching bases.
   - **Strengths**: The script is structured with functions and includes docstrings describing each function and workflow, which improves readability.
   - **Suggested Improvements**: 
     - **Error Handling**: The script ran into a `FileNotFoundError` because `../Results/DNA_seq.txt` was not found. Adding checks to confirm the output directory exists or using `os.makedirs()` would prevent this error.
     - **Comments**: Each alignment function could use additional comments explaining specific steps, especially within `calculate_score` to clarify the purpose of the alignment and scoring.

2. **oaks_debugme.py**:
   - **Purpose**: This script filters a CSV file for oak species (Quercus genus) and writes these to an output CSV.
   - **Strengths**: The script is structured with functions, and it uses docstrings to explain each function’s role, including `is_an_oak` with detailed explanations and test cases.
   - **Suggested Improvements**:
     - **Error Handling**: Consider adding error handling around file reading and writing to manage potential issues with file paths.
     - **Comments**: While the main steps are clear, adding comments explaining the filtering logic and the purpose of each section would improve clarity, especially for functions dealing with CSV operations.

### General Code Suggestions

- **Consistency**: Ensure consistent indentation and spacing across scripts to enhance readability.
- **Error Handling**: Using `try` and `except` blocks for file operations and command-line arguments would make the scripts more robust.
- **Docstrings**: Both scripts include docstrings, but `align_seqs_fasta.py` is missing one for some functions. Consider adding consistent docstrings to clarify functionality.

---