#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: k-Mer Composition
"""

# import statements
from sys import argv

def parse_file(file):
    """ parses an text file document

    file: list; every line of the document as a item
    :return:
    kmer_length: int; the length of the k-mer
    dna: str; DNA sequence
    """
    kmer_length= int(file[0])
    dna=''.join(file[1:]).replace('\n', "")
    return dna, kmer_length

def creating_kmers(dna,kmer_length):
    """Counts the occurrence of the possible k-mers

    dna: str; A DNA sequence
    kmer_length: int; the lenght of the k-mers
    return: list; list of all unique k-mers
    """
    i = 0
    k = kmer_length
    length = len(dna) + 1
    kmers=[]
    while k != length:
        kmer = dna[i:k]
        i += 1
        k += 1
        kmers.append(kmer)
    kmers=list(set(kmers))
    return kmers

def writing_txt_file(kmers):
    """Writing an output file

    kmers: list; returns a list containing all k-mers
    return: None
    """
    with open ('result.txt', 'w') as file:
        for kmer in kmers:
            file.write(str(kmer) + '\n')


def main():
    """This is the main function of the script"""
    file= open(argv[1]).readlines()
    dna,kmer_length=parse_file(file)
    kmers=creating_kmers(dna, kmer_length)
    print(kmers)
   # writing_txt_file(kmers)

if __name__ == "__main__":
    main()