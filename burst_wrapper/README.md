Requires a linux operating system and Python 3.

**Objective:** Complete a convenience wrapper script that runs BURST as a 3rd party software tool for performing alignment, assembly, phylogeny inference, secondary structure prediction, etc. 

First, verified that I could call the BURST command from within python using the following syntax: `./burst -r ref.fna -q query.fna --taxonomy taxonomy.txt -o output.txt`. 

Using [this](https://github.com/seqan/lambda/wiki/BLAST-Output-Formats) format guide for the BLAST output, modified `template.py` into `modified_template.py` to print out the...
- Return value, standard output, standard error of the BURST run. 
- Fraction of the original input query sequences that had a match in the database at 97% or above. 
- Most common bacterial species in the query set. 
- Average percent similarity of the matches, excluding those below 97%. 

Verified that the program runs with `python wrapper.py -q query.fna -r ref.fna -t taxonomy.txt -c ./burst -o output.txt -V`. 

The print out has been logged to `log.txt`.
