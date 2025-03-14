
# Final CMEE Bootcamp Assessment: Tianye

- The project structure followed a logical organization with separate folders for `code`, `data`, `results`, and `sandbox`. 
- Weekly directory and file naming was inconsistent - need more attention to detail on that front.
- `.gitignore` was missing, which led to unnecessary files like `.DS_Store` being tracked in the repository. Adding a `.gitignore` would have improved the workflow and kept the repository cleaner... You might [find this useful](https://www.gitignore.io).
- The `README` files were present in most directories but lacked sufficient detail in earlier weeks. By Week 4, the `README.md` provided comprehensive instructions, usage details, and context for the scripts - a notable improvement. You could have included versions of languages and dependencies/packages used. Also check out [this resource](https://github.com/jehna/readme-best-practices).
- Results folders were appropriately kept empty in most weeks, except for Week 2, where a pre-generated results file was found. This was rectified in subsequent weeks.
- Week 4 introduced clearer documentation and LaTeX integration.

## Week 1
- Shell scripts were generally functional. However, comments in scripts like `CountLines.sh` and `csvtospace.sh` could have been more detailed.

## Week 2
- Python script functions were well-defined, and file paths were correctly referenced.
- Missing dependencies (`fuzzywuzzy`) caused scripts like `oaks_debugme.py` to fail initially. This issue was documented and resolved through the use of virtual environments in Week 3. Avoid using extra packages...
-  You could have formatted the output of certain scripts to be  more neat / organised / informative (compare with my solutions) -- for example `lc1.py` is perfectly functional, but the formatting of the output could have been improved.

## Week 3
- OK. Code could have been more modular in some cases.


## Week 4
-  Florida's temperature autocorrelation analysis (`TAutoCorr.R` and `Florida.R`) stood out for their clarity and use of statistical methods. However, the LaTeX report could have been enhanced with more detailed captions for figures and tables.
- No significant errors, but results could be better commented in the output files for clarity.

- `Florida.R` efficiently calculated correlations and performed permutation tests. The use of base R functions minimized unnecessary dependencies - good. compare with provided solution.
- `TAutoCorr.R` included detailed comments and clear usage instructions. The code  was reasonably efficient , providing a correct answer to the question. The  provided statistical and biological/ecological interpretations in the report could have been stronger; has a somewhat weak conclusion at the end. Consider how the script could have been further optimized for performance with larger datasets.
- The report was well-organized. The figures lacked descriptive captions, and there were minor typographical errors.
- Overall, your Groupwork practicals were all in order, and your group did well in collaborating on it based on the commit/merge/pull history. Check the groupwork feedback pushed to your group repo for more details.

## Git Practices

- The `git log` analysis showed consistent contributions across weeks. While Week 1 had fewer meaningful commits, subsequent weeks included detailed commit messages and incremental improvements.
- The repository size of ~2.35 MiB is reasonable.

## Recommendations for Future Work

1. `.gitignore` - will prevent unnecessary files like `.DS_Store` from being committed.
2. Keep README files consistent in format and ensure they include a summary, setup instructions, and usage details.
3. **Commit Messages**: While commit quality improved, more descriptive messages about changes would be helpful.

---

## Overall Assessment

Well done.

You demonstrated improvement in code quality, documentation, and workflow management throughout the weeks. 

Some of your scripts retained fatal errors which could nave been easily fixed - work on being more vigilant and persistent in chasing down errors the future.

Commenting could be improved -- you are currently erring on the side of overly verbose comments at times (including in your readmes), which is nonetheless better than not commenting at all, or too little! This will improve with experience, as you will begin to get a feel of what is ``common-knowledge'' among programmers, and what stylistic idioms are your own and require explanation. In general though, comments should be written to help explain a coding or syntactical decision to a user (or to your future self re-reading the code!) rather than to describe the meaning of a symbol, argument or function (that should be in the function's docstring in Python for example).

It was a tough set of weeks, but I believe your hard work in them has given you a great start towards further training, a quantitative masters dissertation, and ultimately a career in quantitative biology! 

### (Provisional) Mark

*66*