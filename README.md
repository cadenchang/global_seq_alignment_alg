Usage: python3 align.py < {inputfile}
Input file in fasta format, contains two sequences of DNA to be aligned

Implements global DNA sequence alignment using a simple similarity scoring scheme with a linear gap penalty:
M for each matched pair in the alignment, m for each mismatch, and g for each gap
