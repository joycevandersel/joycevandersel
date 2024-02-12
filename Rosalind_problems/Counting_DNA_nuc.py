#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind Counting DNA Nucleotides
"""

# import statements
from sys import argv

#function
def count_nuc (dna_seq):
    """Return the nucleotide counts of a sequence

    dna_seq: str, DNA sequence
    nuc_count: list with four numbers
    """
    nuc_count= []
    nuc_count.append(dna_seq.count('A'))
    nuc_count.append(dna_seq.count('C'))
    nuc_count.append(dna_seq.count('G'))
    nuc_count.append(dna_seq.count('T'))
    return (nuc_count)

if __name__ == "__main__":
    """this is the main function of the script"""
    dna_seq=open(argv[1]).read()
    nuc_count = count_nuc(dna_seq)
    print(*nuc_count, sep=' ')
