#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Student number: 1091565
Script to make a hierarchical clustering of a csv file. Each row are data
points of a variable (for example gene expression at time points)

The script is able to use different settings, the initial distance can be
calculated using Euclidean or Correlation distance. For calculating the
distance from the made cluster to the other data points/clusters either
separation distance or averge distance can be chosen.

In the main function the files are read and the pipeline of this script is
implemented. All the answers are printed in a separated function.
"""

# import function
import time
from copy import deepcopy
from math import sqrt


def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try:  # will fail on header line in file
            data_points.append(list(map(float, items[1:])))  # skip label
        except ValueError:  # must be the header
            continue
    return data_points


def create_dis_matrix(data_points, dis_func):
    """ Creates a one-sided distance matrix

    :param data_points:list of list; every list is a data point consisting of
            several entries per points.
    :param dis_func: name of function; this can either be Euclidean or
            correlation based.
    :return: list of list one-sided matrix and a list of list of column or row
    names.

    With the dis_func variable this function became universal
    """
    matrix = []
    for i in range(len(data_points)):
        row = []
        for j in range(len(data_points)):
            if i == j:
                break
            elif i > j:
                distance = dis_func(data_points[i], data_points[j])
                row.append(distance)
        matrix.append(row)
    col_names = ['g' + str(name) for name in range(1, len(matrix[-1]) + 2)]
    row_names = ['g' + str(name) for name in range(2, len(matrix) + 1)]
    matrix_names = [row_names, col_names]
    return matrix[1:], matrix_names  # first is empty


def get_euc_dis(pointa, pointb):
    """ Calculating the Euclidean distance between two points (a and b)

    :param pointa: list; containing all values of data point one.
    :param pointb: list; containing all values of data point two.
    :return: int; a distance between a an b.
    """
    points = list(zip(pointa, pointb))
    dis = 0
    for point_tuple in points:
        value_a, value_b = point_tuple
        dis += (value_a - value_b) ** 2
    dis = round(sqrt(dis), 1)
    return dis


def get_corr_dis(a, b):
    """Calculating the correlation distance between two points (a and b)

    :param a:list; containing all values of data point one.
    :param b:list; containing all values of data point two.
    :return: int; a distance between a an b.
    """
    d = len(a)
    top = (d * sum([a[i] * b[i] for i in range(len(a))])) - (sum(a) * sum(b))
    bot = (sqrt(d * sum([num ** 2 for num in a]) - sum(a) ** 2)) * \
          (sqrt(d * sum([num ** 2 for num in b]) - sum(b) ** 2))
    dis = round(1 - top / bot, 5)
    return dis


def create_cluster(new_matrix, matrix_names, new_dis_func, clus):
    """Clusters data points based on the made distance matrix.

    :param new_matrix: list of lists; one-sided distance matrix.
    :param matrix_names: list of two list; first list containing all the row
            names, the second containing the column names.
    :param new_dis_func: name of function; this can either be
    'calculate_new_dis_matrix_sep' or 'calculate_new_dis_matrix_avg'
    :param clus: int; The number of clusters that are going to be made
    :return: list; each item are all gene that have been combined, the length
        of the list is the number of clusters.

    - This function uses "find_smallest_value", "calculate_new_dis_matrix" and
    "update_names_list" until the length of column names is equal to clus.
    - With the 'new_dis_func' it can be specified how to distance gets
    calculated between
    """
    while len(matrix_names[1]) != clus:
        closest = find_smallest_value(new_matrix)
        new_matrix = new_dis_func(closest, new_matrix, matrix_names)
        matrix_names = update_names_list(matrix_names, closest)
        if len(new_matrix[0]) == 0:
            new_matrix = list(filter(None, new_matrix))
            del matrix_names[0][0]
    return matrix_names[1]


def find_smallest_value(matrix):
    """ Finds the coordinates of the smallest value in a one-sided matrix.

    :param matrix: list of lists; one-sided distance matrix.
    :return: list of two ints; One being the row coordinate and one of the
        column coordinate of the location in the matrix [row, column]

    - This function favours the first minimal value it comes across to save
    time.
    - Something to consider with these coordinates is that due to the matrix
    being one sided the row and column index are not interchangeable. Row
    to column is plus one.
    """
    small_row_val = []
    values = []
    for row in matrix:
        small_row_val.append((min(row), row.index(min(row))))
        values.append(min(row))
    coor_small_dis = []
    coor_small_dis.append(small_row_val.index(min(small_row_val)))
    for tup in small_row_val:
        if min(values) in tup:
            coor_small_dis.append(int(tup[1]))
            break
    return coor_small_dis


def calculate_new_dis_matrix_sep(closest, copy_matrix, matrix_names):
    """Creates new row and removes rows and column from the matrix

    :param closest: list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column].
    :param copy_matrix: list of lists; one-sided distance matrix.
    :param matrix_names: list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column].
    :return: list of lists; updated one-sided distance matrix.

    The two points with the closest distance are combined and the new distance
    is calculated based on separation between the combined points and the
    remaining points. This is stored in a list. Then, the rows and columns from
    the closest points (and their complementary row or column) are removed. And
    finally the stored row gets added to the matrix.
    """
    new_row = []
    row_idx, col_idx = closest
    dis_between = copy_matrix[row_idx][col_idx]
    for i in range(len(copy_matrix[-1]) + 1):
        new_dis = 0
        if i != col_idx and i != row_idx + 1:
            try:
                dis1 = copy_matrix[row_idx][i]
            except IndexError:
                dis1 = copy_matrix[i - 1][row_idx + 1]
            try:
                if i - 1 < 0:
                    dis2 = copy_matrix[col_idx - 1][i]
                else:
                    dis2 = copy_matrix[i - 1][col_idx]
            except IndexError:
                dis2 = copy_matrix[col_idx - 1][i]
            new_dis = (dis1 + dis2 - dis_between) / 2
            new_row.append(round(new_dis, 1))
    new_matrix = remove_cluster_colums_rows(closest, copy_matrix)
    new_matrix.append(new_row)
    return new_matrix


def calculate_new_dis_matrix_avg(closest, copy_matrix, matrix_names):
    """Creates new row and removes rows and column from the matrix
    
    :param closest: list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column].
    :param copy_matrix: list of lists; one-sided distance matrix.
    :param matrix_names: list of two ints; One being the row coordinate and one 
        of the column coordinate of the location in the matrix [row, column].
    :return: list of lists; updated one-sided distance matrix.
    
    The two points with the closest distance are combined and the new distance
    is calculated based on average distance between the combined points and the
    remaining points. This is stored in a list. Then, the rows and columns from
    the closest points (and their complementary row or column) are removed. And
    finally the stored row gets added to the matrix.
    """""
    new_row = []
    row_idx, col_idx = closest
    for i in range(len(copy_matrix[-1]) + 1):
        new_dis = 0
        temp = 0
        if i != col_idx and i != row_idx + 1:
            # print(i)
            try:
                temp += copy_matrix[row_idx][i]
            except IndexError:
                temp += copy_matrix[i - 1][row_idx + 1]
            try:
                if i - 1 < 0:
                    temp += copy_matrix[col_idx - 1][i]
                else:
                    temp += copy_matrix[i - 1][col_idx]
            except IndexError:
                temp += copy_matrix[col_idx - 1][i]
            new_clus = matrix_names[1][col_idx].count('g') + \
                       matrix_names[1][row_idx + 1].count('g')
            old_clus = matrix_names[1][i].count('g')
            new_dis = (1 / (old_clus * new_clus)) * temp
            new_row.append(round(new_dis, 1))
    new_matrix = remove_cluster_colums_rows(closest, copy_matrix)
    new_matrix.append(new_row)
    return new_matrix


def remove_cluster_colums_rows(closest, matrix_copy):
    """Removes the distance values from the two closest points from the matrix.

    :param closest:list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column]
    :param matrix_copy: list of lists; one-sided distance matrix.
    :return: copy_matrix: list of lists; of reduced one-sided distance matrix.

    Due to the matrix being one sided the row and column index are not
    interchangeable. Row to column is plus one and vice versa.
    """
    row_idx, col_idx = closest
    # remove the columns
    for rows in matrix_copy:
        # one sided matrix
        if row_idx + 1 < col_idx and len(rows) > col_idx and len(rows) \
                <= row_idx + 1:
            del rows[col_idx]
            del rows[row_idx + 1]
        elif row_idx + 1 < col_idx and len(rows) > row_idx + 1:
            del rows[row_idx + 1]
        elif row_idx + 1 > col_idx and len(rows) > col_idx and len(
                rows) <= row_idx + 1:
            del rows[col_idx]
        elif row_idx + 1 > col_idx and len(rows) > row_idx:
            del rows[col_idx]
            del rows[row_idx]

    # remove the rows
    if col_idx != 0:  # one-sides matrix so col_idx does not exist for col = 0
        if row_idx > col_idx - 1:
            del matrix_copy[row_idx]
            del matrix_copy[col_idx - 1]
        else:
            del matrix_copy[col_idx - 1]
            del matrix_copy[row_idx]
    else:
        del matrix_copy[row_idx]
    return matrix_copy


def update_names_list(matrix_names, closest):
    """ Removes the names of the two closest points and adds these combined.

    :param matrix_names: list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column]
    :param closest: list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column]
    :return: list of two ints; One being the row coordinate and one of
        the column coordinate of the location in the matrix [row, column]

    This function removes the names of the two points from the row and column
    list, and adds these in the following format back to the row and column.
    'name1, name2'. This makes it easy to print.
    """
    row_idx, col_idx = closest
    row_names, col_names = matrix_names
    gene1 = col_names[col_idx]
    gene2 = row_names[row_idx]
    if gene1 in row_names:
        row_names.remove(gene1)
    if gene2 in row_names:
        row_names.remove(gene2)
    if gene1 in col_names:
        col_names.remove(gene1)
    if gene2 in col_names:
        col_names.remove(gene2)
    row_names.append(str(gene1 + ", " + gene2))
    col_names.append(str(gene1 + ", " + gene2))
    return matrix_names


def pipeline(file, dis_func, new_dis_func, clus):
    """

    :param file: list of lines; an file opened using read lines.
    :param dis_func: name; the function name of the methode to calculate the
            distance.
    :param clus: int; The number of clusters that are going to be made
    :return: list of list one-sided matrix and a list of cluster names.
    """
    point = csv_parser(file)
    #point = transpose(point)
    dis_matrix, matrix_names = create_dis_matrix(point, dis_func)
    matrix_names = create_cluster(deepcopy(dis_matrix), matrix_names,
                                  new_dis_func, clus)
    return dis_matrix, matrix_names


def print_question_answers(Q1_dis_matrix, Q1_clusters, Q2a_dis_matrix,
                           Q2a_clusters, Q2b_dis_matrix, Q2b_clusters,
                           Q7_clusters, Q8_clusters):
    """ This function prints the answers to the questions

    :param Q1_dis_matrix: list of list one-sided matrix of 'jp_fig10_1a.csv'
        using Euclidean distance.
    :param Q1_clusters: list of cluster names of 'jp_fig10_1a.csv' using
        Euclidean distance.
    :param Q2a_dis_matrix: list of list one-sided matrix of
        'example_gene_expression_four_genes.csv' using Euclidean distance.
    :param Q2a_clusters: list of cluster names of
        'example_gene_expression_four_genes.csv' using Euclidean distance.
    :param Q2b_dis_matrix: list of list one-sided matrix of
        'example_gene_expression_four_genes.csv' using Correlation distance.
    :param Q2b_clusters: list of cluster names of
        'example_gene_expression_four_genes.csv' using Correlation distance.
    :param Q7_clusters: list of cluster names of 'proteomics_data.csv' using
        Correlation distance.
    :param Q8_clusters:
    :return: None
    """
    print("Question 1a:")
    for row in Q1_dis_matrix:
        print('\t', row)
    print("Question 1b:")
    for i, clus in enumerate(Q1_clusters):
        print('\t cluster {}: {}'.format(i + 1, clus))
    print("Question 2a:")
    for row in Q2a_dis_matrix:
        print('\t', row)
    print("Question 2b:")
    for row in Q2b_dis_matrix:
        print('\t', row)
    print("Question 2c:")
    for i, clus in enumerate(Q2a_clusters):
        print('\t cluster {}: {}'.format(i + 1, clus))
    print("Question 2d:")
    for i, clus in enumerate(Q2b_clusters):
        print('\t cluster {}: {}'.format(i + 1, clus))
    print("Question 7:")
    for i, clus in enumerate(Q7_clusters):
        print('\t cluster {}: {} genes'.format(i + 1, len(clus)))
    print("Question 8:")
    for i, clus in enumerate(Q1_clusters):
        print('\t cluster {}: {}'.format(i + 1, clus))


def main():
    """This is the main function of the script"""
    start_time = time.time()
    # Question 1
    file = open(
        '../../../PycharmProjects/opdracht5/jp_fig10_1a.csv').readlines()
    Q1_dis_matrix, Q1_clusters = pipeline(file, get_euc_dis,
                                          calculate_new_dis_matrix_sep, clus=3)

    # Question 2
    file = open(
        '../../../PycharmProjects/opdracht5/example_gene_expression_four_genes.csv').readlines()
    Q2a_dis_matrix, Q2a_clusters = pipeline(file, get_euc_dis,
                                            calculate_new_dis_matrix_sep,
                                            clus=2)
    # Needs to be recalled due to shadowing
    file = open(
        '../../../PycharmProjects/opdracht5/example_gene_expression_four_genes.csv').readlines()
    Q2b_dis_matrix, Q2b_clusters = pipeline(file, get_corr_dis,
                                            calculate_new_dis_matrix_sep,
                                            clus=2)

    # Question 7
    file = open(
        '../../../PycharmProjects/opdracht5/proteomics_data.csv').readlines()
    Q7_dis_matrix, Q7_clusters = pipeline(file, get_corr_dis,
                                          calculate_new_dis_matrix_sep, clus=3)

    print("Question 7:")
    for i, clus in enumerate(Q7_clusters):
        print('\t cluster {}: {} genes'.format(i + 1, clus))
    #Question 8
    file = open(
        '../../../PycharmProjects/opdracht5/jp_fig10_1a.csv').readlines()
    Q8_dis_matrix, Q8_clusters = pipeline(file, get_euc_dis,
                                          calculate_new_dis_matrix_avg, clus=3)

    # printing
    print_question_answers(Q1_dis_matrix, Q1_clusters, Q2a_dis_matrix,
                           Q2a_clusters, Q2b_dis_matrix, Q2b_clusters,
                           Q7_clusters, Q8_clusters)
    end_time = time.time()
    print(end_time - start_time)

def transpose(datapoints):
    a = list(zip(*datapoints))
    z = [list(x) for x in a]
    return z


if __name__ == "__main__":
    main()
