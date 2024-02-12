#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Script to solve Rosalind challenge: Breadth-First Search
"""

# import statements
from sys import argv

def parse_file(file):
    """ Parsing an edge list format
    file:list; containing every line of a edge list format
    return:
        nodes:int; The number of points
        connections:list; direction of the connections between the points

    In this it is important to recognise what the data looks like. This is a directed graph,
    with the nodes labeled from 1 to n.
    """
    nodes_edges=file[0].replace('\n', "").split(' ')
    nodes=int(nodes_edges[0])
    con=file[1:]
    connections='\n'.join(con).replace(' ',"-").split('\n')
    return(nodes, connections)

def creat_default_dic(nodes):
    """Here a dictionary is build for each of the nodes to get a connection value

    nodes:int; The number of points
    returns: dict; the key represent a node, the value is set to zero"""
    distance_dic={}
    for x in range(1,nodes+1):
        distance_dic[x] = 0
    return(distance_dic)


def going_true_nodes(con,dis_dic,nodes,start):
    """ Calculating the distance of nodes from a starting point

    :input:
    con:list; direction of the connections between the points.
    dis_dic:dict; the key represent a node, the value is set to zero.
    nodes:int; The number of points.
    start:list; which node the function should start in a list.

    :return:
    dis_dic:dict; the key represents a node, the value is the distance from start to node.
    (-1 represents no options).
    """
    list_of_keys=[start]
    run=nodes+1
    for x in range(1,run):
        list_of_keys, dis_dic = distance_calculator(con, dis_dic, list_of_keys, x)
    list_of_keys = [key for (key, value) in dis_dic.items() if value == 0]
    for key in list_of_keys:
        dis_dic[key] = -1
        dis_dic[start]=0
    return(dis_dic)


def distance_calculator(con,dis_dic,list_of_keys, times):
    """ Calculating the distance of nodes from a starting point

    :input:
    con:list; direction of the connections between the points.
    dis_dic:dict; the key represent a node, the value is set to zero.
    list_of_keys:list; list of all possible connecting nodes
    times:int; how many previous connections were needed to get to this point

    :return:
    list_of_keys:list; list of all possible connecting nodes
    dis_dic:dict; the key represents a node, the value is the distance from start to node.
    (-1 represents no options).
    """

    for i in list_of_keys:
        for node in con:
            coor = node.split('-')
            a = int(coor[0])
            b = int(coor[1])
            if i == a and dis_dic[b] == 0:
                 dis_dic[b] += times
    list_of_keys = [key for (key, value) in dis_dic.items() if value == times]
    return (list_of_keys, dis_dic)


def getting_distace(dis_dic):
    """creates a list of the distance

    dic_scores: dict; keys are the nodes, values are the distance from start node to node
    return: list; returns a list containing all values of a dictionary
    """
    distance=list(dis_dic.values())
    return (distance)

def writing_txt_file(distance):
    """Writing an output file

    kmers: list; list containing all distances
    return: None
    """
    with open ('result.txt', 'w') as file:
        for dis in distance:
            file.write(str(dis) + ' ')


def main():
    """this is the main function of the script"""
    file= open(argv[1]).read().strip().split('\n')
    start=1
    nodes, connections = parse_file(file)
    dis_dic=creat_default_dic(nodes)
    going_true_nodes(connections,dis_dic, nodes,start)
    dis=getting_distace(dis_dic)
    writing_txt_file(dis)


if __name__ == "__main__":
    main()