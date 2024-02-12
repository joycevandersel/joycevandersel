#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: Constructing a De Bruijn Graph
"""

# import statements
from sys import argv

def making_union(kmers):
    """Finding all edges between the kmers
    :input:
        kmers: list; every item is a DNA price called S
    return:
        kmer_Set: list; list containing tuples of (k,k)
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    kmer_Set = set()
    for kmer in kmers:
        rev_kmer="".join(complement.get(nuc, nuc) for nuc in reversed(kmer))
        kmer_Set.add(kmer)
        kmer_Set.add(rev_kmer)
    kmer_Set=sorted(kmer_Set)
    return (kmer_Set)

def retrieving_edges(kmer_set):
    """Finding all edges between the kmers
    :input:
        kmer_set: set; a set containing all S and S^rc
    return:
        edges: list; list containing tuples of (k,k)
    """
    list1=[]
    list2=[]
    for kmer in kmer_set:
        x=len(kmer)-1
        list1.append(kmer[0:x])
        x=x+1
        list2.append(kmer[1:x])
    edges=list(zip(list1,list2))
    return(edges)

def writing_txt_file(edges):
    """Writing an output file

    edges: list; list containing tuples of (s,s^rc)
    return: None
    """
    #for edge in edges:
    with open ('result.txt', 'w') as file:
        for edge in edges:
            if len(edge[0]) != 0:
                file.write("({}, {})\n".format(edge[0], edge[1]))


def main():
    """this is the main function of the script"""
    kmers = open(argv[1]).read().split("\n")
    kmers_set=making_union(kmers)
    edges=retrieving_edges(kmers_set)
    writing_txt_file(edges)

if __name__ == "__main__":
    main()
