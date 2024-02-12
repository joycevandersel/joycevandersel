#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: Introduction to Pattern Matching
"""

# import statements
from sys import argv
import time


def extending_dic(dna_string,pathway_dic, counter):
    """ making an adjacency list based on pattern matching

    :input:
        dna_string: str; a string containing DNA nucleotides
        pathway_dic:dic; a dictionary containing the excisting connections for
            node each
        counter:int; keeps track of how many nodes exist

    :return:
        pathway_dic:dic; a dictionary containing the excisting connections for
            each node
        counter:int; keeps track of how many nodes exist

    This function is called multiple times each with a different DNA string.
    The counter and dictionary are returned to this function with a new string.
    Every time the dictionary get extended with new connections from starting
    node 1. Being DNA there can be no more than four connections in any value
    """
    node=1
    for x in range(len(dna_string)):
        if any(dna_string[x] in sl for sl in pathway_dic[node]):
            for item in pathway_dic[node]:
                if dna_string[x] in item:
                    node = item[0]
        else:
            pathway_dic[node].append([counter, dna_string[x]])
            node = counter
            pathway_dic[counter] = []
            counter = counter + 1
    return (pathway_dic, counter)

def cleaning_for_end_nodes(pathway_dic):
    """ cleaning the dictonary of empty lists
    :input:
    pathway_dic: dict; key is the start and end node, value are nested list
    containing connected node followed by the name of the branch

    :return:
    pathway_dic: dict; key is the start node, value are nested list containing
    connected node followed by the name of the branch
    """
    empty_nodes=[]
    [empty_nodes.append(k) for k, v in pathway_dic.items() if v == []]
    [pathway_dic.pop(i) for i in empty_nodes]
    return (pathway_dic)

def printing_results(pathway_dic):
    """printing the results

    :input:
    pathway_dic: dict; key is the start node, value are nested list containing
    connected node followed by the name of the branch
    :return: None
    """
    for k, v in pathway_dic.items():
        for item in v:
            print(k, *item, sep=' ')



def main():
    """this is the main function of the script"""
    start_time=time.time()
    dna_list = open(argv[1]).read().split('\n')
    dna_list = sorted(dna_list, key=len, reverse=True)
    pathway_dic={1:[]}
    counter=2
    for item in dna_list:
        pathway_dic, counter=extending_dic(item,pathway_dic, counter)
    pathway_dic=cleaning_for_end_nodes(pathway_dic)
    printing_results(pathway_dic)
    end_time = time.time()
    print(end_time-start_time)

if __name__ == "__main__":
    main()