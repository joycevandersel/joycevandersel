#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: Connected Components
"""

# import statements
from sys import argv


def parse_file(file):
    """ Parsing an edge list format
    file:list; containing every line of a edge list format
    return:
        nodes:int; The number of points
        connections:list; direction of the connections between the points

    In this it is important to recognise what the data looks like. This is a
    directed graph, with the nodes labeled from 1 to n.
    """
    nodes_edges = file[0].replace('\n', "").split(' ')
    nodes = int(nodes_edges[0])
    con = file[1:]
    connections = '\n'.join(con).replace(' ', "-").split('\n')
    connections = [x for x in connections if x.strip()]
    return (nodes, connections)


def creating_dict(con, nodes):
    """Makes a dictonary containig all links

    :input:
        con:list; direction of the connections between the points ['int-int']
        nodes:int; The number of points

    :return:
        nodes_dic:dict; the key represents a node, the value is a list of nodes
                it is connected to. Here every combination is present. meaning
                a conection between 1 and 2 is also between 2 and 1 in the
                disctionary.
    """
    nodes_dic={}
    for x in range(1,nodes+1):
        nodes_dic[x] = [x]
    for node in con:
        coor = node.split('-')
        a = int(coor[0])
        b = int(coor[1])
        for x in range (1,nodes+1):
            if x==a:
                nodes_dic[x].append(b)
            if x==b:
                nodes_dic[x].append(a)
    return(nodes_dic)

def conecting_dict(dic):
    """Combines all values from dictionary into unique sets

    :input:
        dic:dict; the key represents a node, the value is a list of nodes
                it is connected to.
    :return:
     out:list of sets; it is a list containing sets of numbers

    This function takes every the first list of a nested list in a while loop
    until it is empty. Then it makes sets of the first item in a list and
    compares this to the already excisting sets. The excisting set gets
    extended or if it does not excist it is added.
    """
    con=list(dic.values())
    graphs = []
    while len(con) > 0:
        first, *rest = con
        first = set(first)
        x = -1
        while len(first) > x:
            x = len(first)
            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2
        graphs.append(first)
        con = rest
    return(graphs)


def main():
    """this is the main function of the script"""
    file = open(argv[1]).readlines()
    nodes, connections = parse_file(file)
    nodes_dic=creating_dict(connections, nodes)
    graphs=conecting_dict(nodes_dic)
    print(len(graphs))


if __name__ == "__main__":
    main()
