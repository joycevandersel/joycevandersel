#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: Generate the k-mer Composition of a String
"""

# import statements
from sys import argv


def parse_file(file):
    """Returns only the nucleotide sequence without header

    dna_seq:list; text-file
    nuc_count: str; dna sequence

    The input file is a text file containing a header followed by a
    dna sequence in a list format.
    """
    dna_seq = ''.join(file[1:]).replace('\n', "")
    return (dna_seq)


def creating_kmers_list(kmer, kmer_length, next_list):
    """Creatin al list of all possible DNA combinations

    kmer: str; an empty string
    kmer_length: int; base length of the k-mer
    next_list: list; an empty list for storing the k-mers
    return: this breaks the loop created and the next_list gets returns

    This is recursive function.
    """
    nuc = ['A', 'T', 'G', 'C']
    if kmer_length==0:
        next_list.append(kmer)
        return
    for i in range(len(nuc)):
        new_kmer= kmer + nuc[i]
        creating_kmers_list(new_kmer,kmer_length-1, next_list)

def creating_dic(kmers_list):
    """creates a dictionary based on a list

    kmers_list: List; of all possible k-mers
    return: dict; keys are the k-mers, values are 0
    """
    sorted_list=sorted(kmers_list)
    kmers_dict={}
    for k in sorted_list:
        kmers_dict[k]=0
    return kmers_dict

def kmers_composition(kmers_dict, dna):
    """Counts the occurrence of the possible k-mers

    kmers_dict: dict; keys are the k-mers, values are 0
    dna: str; A DNA sequence
    return: dict; keys are the k-mers, values are the occurrence
    """
    i = 0
    k = 4
    length = len(dna) + 1
    while k != length:
        kmer = dna[i:k]
        i += 1
        k += 1
        kmers_dict[kmer] += 1
    return(kmers_dict)

def getting_scores(dic_scores):
    """creates a list of the occurrence score

    dic_scores: dict; keys are the k-mers, values are the occurrence
    return: list; returns a list containing all values of a dictionary
    """
    scores=list(dic_scores.values())
    return (scores)

def writing_txt_file(scores):
    """Writing an output file

    scores: list; returns a list containing all scores
    return: None
    """
    print(scores)
    with open ('result.txt', 'w') as file:
        ' '.join(map(str, scores))
        #for score in scores:
            #file.write(str(score) + ' ')

def main():
    """this is the main function of the script"""
    file = open(argv[1]).readlines()
    kmer_length = 4
    dna_seq = parse_file(file)
    kmer_list = []
    creating_kmers_list('',kmer_length, kmer_list)
    kmers_dict=creating_dic(kmer_list)
    dic_score=kmers_composition(kmers_dict, dna_seq)
    scores=getting_scores(dic_score)
    writing_txt_file(scores)


if __name__ == "__main__":
    main()



