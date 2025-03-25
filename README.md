# minERVa
MinERVa is the provisory name for a program taht runs on *nix systems and is meant to automatize analysis of the complete diversity of a single species of endogenous retrovirus (ERV) contained in its host genome. An enormous amount of long-read genomic data is being made available every day, and up to 10% of host's genomes are composed of ERVs. Cannonically, ERVs have been known by many names: "junk DNA" was the first, and nowadays "genomic fossil" is the most used term. These terms indicate a trend in ERV thought: that they are inactive, remnants of past organisms. It is considered that ERVs are slow-evolving entities, affected only by genetic drift. My greatest aim is to aid researchers worldwide in proving otherwise.
This program aims to facilitate the analysis of all copies of a single ERV and point to possible roles in host biology. So far, the output of the automatized script is the region surrounding the ERV, which can be applied to an annotation program in order to define the beggining and end of the ERV itself. 
The accompanying document is a manual that can be used to guide further studies, and my aim with new versions is to automatize everything that is in it and more. My end goal is to determine a single ERVs diversity, organizing the data in a manner that is more easily interpretable by researchers, and didatic to those that are beginning in ERV studies. Diverse analysis will be performed in order to build this report and can be done independently and their outputs can be extracted for future use: annotation of copies' contents, proportion of copy types, datings of insertions, open-reading frame (ORF) and stop codon analysis, and even phylogenetic analysis. The objetctive is to have an easy overview of any ERV of interest, allowing for better targeting for future studies.

Installation

Usage
Inputs and Outputs
-report - entrega um arquivo PDF descrevendo resultados visualmente 
-phylogeny - entrega um arquivo .treefile
-alignment - entrega um arquivo .fasta alinhado e trimado
-notrimalignment - entrega o arquivo .fasta sem trimar
-karyotype entrega - o arquivo do ChromoMap
-annotation - entrega uma tabela com as anotações (buscar um programa de linha de comando para anotação).
-orf - entrega uma tabela com as ORFs encontradas
-stopcodons - entrega uma tabela com o número de stop codons em cada frame de leitura de cada cópia

Running minERVa

Notes

Citation
A scientific publication fully describing this pipeline is being prepared. Meanwhile, feel free to cite this GitHub repo. Primary references for used dependecies should also be cited:
