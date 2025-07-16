# MINERVA
# Introduction
MinERVa is a program that runs on *nix systems meant to automatize analysis of the complete diversity of a single species of endogenous retrovirus (ERV) contained in its host genome. An enormous amount of long-read genomic data is being made available every day, and up to 10% of host's genomes are composed of ERVs. Cannonically, ERVs have been known by many names: "junk DNA" was the first, and nowadays "genomic fossil" is the most used term. These terms indicate a trend in ERV thought: that they are inactive, remnants of past organisms. It is considered that ERVs are slow-evolving entities, affected only by genetic drift. My greatest aim is to aid researchers worldwide in proving otherwise.
This program aims to facilitate the analysis of all copies of a single ERV and point to possible roles in host biology. So far, the output of the automatized script is the region surrounding the ERV, which can be applied to an annotation program in order to define the beggining and end of the ERV itself. 
The accompanying document is a manual that can be used to guide further studies, and my aim with new versions is to automatize everything that is in it and more. My end goal is to determine a single ERVs diversity, organizing the data in a manner that is more easily interpretable by researchers, and didatic to those that are beginning in ERV studies. Diverse analysis will be performed in order to build this report and can be done independently and their outputs can be extracted for future use: annotation of copies' contents, proportion of copy types, datings of insertions, open-reading frame (ORF) and stop codon analysis, and even phylogenetic analysis. The objetctive is to have an easy overview of any ERV of interest, allowing for better targeting for future studies.

# Installation
The minERVa program requires no installation. However, it requires the certain programs to be available in the directory. Following are the programs required and a link to facilitate installation.

1. BLAST+: https://www.ncbi.nlm.nih.gov/books/NBK569861/
2. MAFFT: https://mafft.cbrc.jp/alignment/software/
3. IQ-TREE3: https://iqtree.github.io/
4. FIGTREE: https://github.com/rambaut/figtree/releases
5. Python packages: os, argparse, csv, re, datetime

Reminder: the programs must be available in the same working directory as minERVa, either add said programs to PATH or use a package manager such as Conda.

# Usage
In order to use minERVa, the user must simply call upon the program in a directory that contains "minERVa_021.py", "retroviridae_references.fasta", the query and genome files. The basic command to activate minERVa is:

python3 minERVa_021.py -genome 'genome_file' -query 'query_file'

This command mines an ERV family according to its LTR ('-query') in one or more genomes (-genome). The outputs at the moment are:
1. the BLAST results related to this run;
2. the location of each copy after filtering;
3. .fasta files containing the sequence of each copy including a genomic context and trimmed by their genomic signature (or inverted repeat) and aligned with a standard retroviral alignment;
4. a ML phylogenetic analysis performed by IQ-TREE, which is automatically opened in FIGTREE for manual evaluation.

Besides the mandatory arguments, some flexibility has already been included in minERVa. The following arguments allow the full use of minERVa:
  
  -h, --help
  show a help message and exit
    
  -genome GENOME [GENOME ...], -g GENOME [GENOME ...], -subject GENOME [GENOME ...]
  path to the input file(s) containing the genomes to be mined in .FASTA format. It is possible to input multiple genomes as a single concatenated file or in different files by specifying each path separated by a space.
  
  -dbready, -ready      
  signals that the file provided by the genome argument is already a database for the BLAST+ program
  
  -query QUERY [QUERY ...], -q QUERY [QUERY ...]
  path to the input file containing the LTR sequence to be used as query in BLAST searches.
  
  -hspcoverage HSPCOVERAGE, -coverage HSPCOVERAGE
  a numerical value between 0 and 100 that specifies the threshold for High-scoring pairs coverage to be used in the BLAST search. Default is 50
  
  -identity IDENTITY, -id IDENTITY
  a numerical value between 0 and 100 that specifies the threshold for identity to be used in the BLAST search. Default is that this threshold is not used       
  
  -context_range CONTEXT_RANGE, -contextr CONTEXT_RANGE
  a numerical value representing the number of base pairs around the hit to remove duplicate results. Default is 20000 base pairs.
  
  -dbkeep, -keep
  signals that the database prepared by this program should not be deleted by the end of the program run. This argument is ignored if the option -dbready is active.
  
  -LTR_range LTR_RANGE, -LTRr LTR_RANGE
  a numerical value representing the number of base pairs around the hit to ensure a perfect trim of the ERV limits. Default is 500 base pairs.
  
  -prefix PREFIX, -run PREFIX, -nick PREFIX
  a string that contains a prefix or nickname for a single minERVa run.

# FUTURE DIRECTIONS:
# Inputs and Outputs
I intend to automatize and flexibilize the analysis of ERVs. In order to do that, I'll consider more input file options, e.g. if the user already has the ERV files and wants only a specific analysis; and more output files, e.g. specific analyses, reports or even images. This is an incomplete list of intended implementations:

-ervsequence - receives a .fasta file containing a set of ERVs as input

-report - outputs a .pdf file with a visual and tabular representation of results

-alignment - outputs a .fasta file of aligned ERVs against a standard retroviral alignment.

-phylogeny - outputs a phylogenetic inference of the mined copies against a standard retroviral alignment.

-karyotype - outputs a figure with the position of the ERVs in relation to the host chromosomes.

-annotation - outputs a table with gene and LTR annotations.

-orf - outputs a table of ORF regions found in the ERVs

-stopcodons - outputs a table containing the stop codon count in each reading frame of each copy

# Notes

# Citation
A scientific publication fully describing this pipeline is being prepared. Meanwhile, feel free to cite this GitHub repo. Primary references for used dependecies should also be cited:
