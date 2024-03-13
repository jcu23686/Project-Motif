#!/bin/bash

# Compile and run the C++ program
g++ gibbs16.cpp -o gb
./gb E.coliRpoN-sequences-16-100nt.fasta 22

# Run the first Python program
python3 PSSM_MOTIF_FINDER.py FruR.txt ecoK12-MG1655.fasta 18

# Run the second Python program
python3 logos.py