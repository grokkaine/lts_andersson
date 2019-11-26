# CRISPR analysis in Baltic metagenomes

This repository contains Andersson project source files, with Sergiu 
Netotea as main and Mikael Huss as co-contributor. Most work is 
summarized in the form of multiple Jupyter notebooks,
in the /notebooks folder, and only some of the code was sourced.

The code was used either from the processing pipeline or from the
notebooks, therefore there is no main script, or command line options,
or even clear inputs and outputs. Instead individual functions are 
called when needed via the notebooks. It can be viewed as a function 
library, but much of is not general enough to qualify, it is for the 
most part internal project code.


Here are the highlights, in no particular order:

- notebooks/index.ipynb - Incomplete index of the notebooks, but it should
be the starting point to finding the most relevant ones. They are grouped by 
processing and analysis. By data analysis should be understood the 
processing of the CRISPR results obtained via CRASS (anything downstream 
from the CRASS program call)


- data processing pipeline, at /pipeline, using snakemake. Most of 
the processing is documented in the processing logs (check the index).


- src/crass_crispr_parser.py - Main file of functions used to parse the 
XML output generated by CRASS. The way to use it is documented in the
"analysis_py" notebooks (check the notebook index).


- src/graph_retrieval.py - Used to generate an older interactive
spacer graph plot, details in notebooks/logs/spacer_graphs.ipynb


- src/pipeline.py and preprocessing.py - Called from the data processing 
pipeline.