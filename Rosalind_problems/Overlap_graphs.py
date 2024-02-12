#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: Overlap Graphs
"""

# import statements
from sys import argv

def creating_list_and_dic (fasta_list):
    """ Creating a header, nucleotide list and nucleotide header dictionary

    :input:
        fasta_list: list; Every item is: header\nsequence
    :return:
        edges_list:nested list; the sublists containing the header
            followed by the last three nucleotides.
        start_dic:dictonary; the key are the last three nucleotides, the
        values are header names that have the key as last three nucleotides.
    """
    edges_list=[]
    start_dic={}
    for item in fasta_list:
        item=item.split('\n')
        seq=''.join(item[1:])
        if seq[0:3] not in start_dic:
            start_dic[seq[0:3]] = [item[0]]
        else:
            start_dic[seq[0:3]].append(item[0])
        edges_list.append([item[0], seq[-3:]])
    return (edges_list, start_dic)


def finding_con(edges_list, start_dic):
    """Makes a set of all possible combinations

    :input:
        edges_list:nested list; the sublists containing the header
            followed by the last three nucleotides.
        start_dic:dictonary; the key are the last three nucleotides, the
        values are header names that have the key as last three nucleotides.
    :return:
        adjecent: set; the set contains tuples consisting of the fasta header
        names (fasta name 1, fasta name 2)

    This function calls 'find links'
    """
    adjecent=set()
    for item in edges_list:
        links=find_links(item, start_dic)
        for i in links:
            if i[0] != i[1]:
                adjecent.add(tuple(i))
    return (adjecent)

def find_links(item, start_dic):
    """
    :input:
        item:list; containing the header followed by the last three nucleotides
        start_dic: dictonary; the key are the last three nucleotides, the
        values are header names that have the key as last three nucleotides.
    :return:
        links:nested list; the sublist contain the two headers that have a
        connection (the 3 nucleotides)
    """
    end = item[1]
    links=[[]]
    if end in start_dic:
        if len(start_dic[end]) == 1:
            if item[0] != start_dic[end]:
                links[-1].append(item[0])
                links[-1].extend(start_dic[end])
                links.append([])
        else:
            l=start_dic[end]
            for x in range(len(l)):
                if item[0] != l[x]:
                    links[-1].append(item[0])
                    links[-1].extend([l[x]])
                    links.append([])
    #the last one is empty
    return(links[:-1])



def writing_txt_file(adjecent):
    """Writing an output file

    adjecent: set; set of tuples (fasta name 1, fasta name 2)
    return: None
    """
    with open ('result.txt', 'w') as file:
        for ad in adjecent:
            if len(ad[0]) != 0:
                file.write("{} {}\n".format(ad[0], ad[1]))

def main():
    """this is the main function of the script"""
    fasta_list = open(argv[1]).read().split('>')
    edges_list, start_dic= creating_list_and_dic(fasta_list[1:])
    adjecent=finding_con(edges_list, start_dic)
    writing_txt_file(adjecent)

if __name__ == "__main__":
    main()