#!/usr/bin/env python3

"""
Author: Joyce van der Sel

Description: this is a script to make a global alignment

Usage: python [script.py]
In the scoring matrix lines containing extra information should contain a # at
the beginning of the line.
"""
#import statement(s)
from sys import argv


def blosum_parser(blosum):
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

# NM: PEP 8 format requires two blank lines between functions, was one
# NM: normally it's advised to have a verb in the function's name
def initial_matrix(seq1, seq2, end_gap_pen):
    """ Creating matrix with end gap set

    seq1: str; first sequence
    seq2: str; second sequence
    end_gap_pen: int; end gap penalty

    :return:
        ini_matrix:list of list; a matrix consiting of a list of list
            containing all end gap penalty

    # NM: Not completely clear
    seq 1 is on are the rows, seq 2 the columns. can be accesed by:
    matrix[row][column]
    """
    ini_matrix = [[0 for colum in range(len(seq2)+1)] \
                  for row in range(len(seq1)+1)]
    for row_idx in range(1, len(seq1)+1):
        # NM: PEP 8 warning
        ini_matrix[row_idx][0] = ini_matrix[row_idx -1][0]+end_gap_pen
    for col_idx in range(1, len(seq2)+1):
        # NM: PEP 8 warning
        ini_matrix[0][col_idx] = ini_matrix[0][col_idx -1]+end_gap_pen
    # NM: parenthesis are not needed here
    return (ini_matrix)


# NM: PEP8 warning about spaces between parameters
# NM: There are more PEP8 warning in other places of the code, I won't mention
# all of them. I wrote some recommendations in my Teams message.
def fil_matrix(seq1, seq2,matrix,gap_pen,order,blosum_matrix):
    """ Filling initial matrix with values

    seq1: str; first sequence.
    seq2: str; second sequence.
    matrix:list of list; a matrix consiting of a list of list
            containing all end gap penalty.
    gap_pen: int; gap penalty.
    order:dict of {res: idx_in_matrix}.
    blosum_matrix: list of lists with similarity scores.

    :return:
        matrix: list of list; a matrix consisting of a list of list containing
            the best possible alignment score per position.
        position:list; The final end position this is where the traceback needs
            to start [row index, colom index].
        traceback: list of list; matrix with each point containing where it
            originated from.

    # NM: I would say that such detailed description is redundant.
    # NM: The code should be self-explanatory. And your code is quite clear,
    # NM: so I would remove this part.
    -This function goes through a matrix row by row (skipping the first value,
    because this was filled all ready). it looks at the previous value
    horizontal, vertical and diagonal up. Applies the appropriate score/penalty
    chooses the best option and set this in the matrix.
    -Next in the tracback matrix the value is the chosen direction. This is
    useful when doing a trace back.
    """
    que=[[row,col] for row in range(1,len(seq1)+1) for col in\
         range(1,len(seq2)+1)]
    traceback = [['' for colum in range(len(seq2)+1)] \
                  for row in range(len(seq1)+1)]
    for position in que:
        optional_values={'dia':0, 'hor':0, 'ver':0}
        ver=position.copy()
        hor=position.copy()
        #Vertical move
        ver[0]=ver[0]-1
        optional_values["ver"]=matrix[ver[0]][ver[1]]+ gap_pen
        #Horizontal move
        hor[1]=hor[1]-1
        optional_values["hor"]=matrix[hor[0]][hor[1]] + gap_pen
        #Diagonal move
        dia=[i-1 for i in position]
        #calls function getting function
        value= getting_value(order, blosum_matrix, dia, seq1,seq2)
        optional_values["dia"]=matrix[dia[0]][dia[1]] + value
        #optional_values.append(matrix[dia[0]][dia[1]] + value)
        best_direction = max(optional_values, key=optional_values.get)
        matrix[position[0]][position[1]]=optional_values[best_direction]
        traceback[position[0]][position[1]]=best_direction
    return (matrix, position, traceback)

# NM: I would be more specific in name, so it's clear which value is being got.
# NM: Something like get_alignment_score.
def getting_value(order, blosum_matrix, dia, seq1, seq2):
    """Looking up value in Blosum matrix

    order:dict of {res: idx_in_matrix}.
    blosum_matrix: list of lists with similarity scores.
    # NM: why is the list of coordinates is called dia?
    dia: list; coordinates to lookup in the blossum matrix [int, int].
    seq1: str; first sequence.
    seq2: str; second sequence.

    :return:
        value: int; alignment score. The closer the AA the higher the score
    """
    AA1 = seq1[dia[0]]
    AA2 = seq2[dia[1]]
    pos = order[AA1]
    pos2 = order[AA2]
    value = blosum_matrix[pos][pos2]
    return(value)


def completing_traceback(traceback, seq1, seq2):
    """Filling first row and colum of the traceback matrix

    traceback:
    seq1: str; first sequence.
    seq2: str; second sequence.
    :return:
        traceback: list of list; matrix with each point containing where it
            originated from.

    the reason these need to be filled is to ensure return to point [0,0]
    """
    for row_idx in range(1, len(seq1)+1):
        traceback[row_idx][0] = 'ver'
    for col_idx in range(1, len(seq2)+1):
        traceback[0][col_idx] = 'hor'
    return traceback

def traceback_matrix(traceback, start, seq1, seq2):
    """ Tracing back trough the matrix to find optimal alignmet

    traceback: list of list; matrix with each point containing where it
            originated from.
    start: list; The final end position this is where the traceback needs
        to start [row index, colom index].
    seq1: str; first sequence.
    seq2: str; second sequence.

    :return:
        align1: list; either an AA or a - meaning a gap the sequence
        matches: list; either a space or a | to indicate two AA machting up
        align2:list; either an AA or a - meaning a gap the sequence
    """
    alignment = [[], [], []]
    while start != [0,0]:
        dir=traceback[start[0]][start[1]]
        row_idx=start[0]-1
        col_idx=start[1]-1
        if dir== "dia":
            start = [i - 1 for i in start]
            alignment[0].append(seq1[row_idx])
            alignment[2].append(seq2[col_idx])
            if seq2[col_idx] == seq1[row_idx]:
                alignment[1].append('|')
            else:
                alignment[1].append('')
        elif dir == "ver":
            start[0]= start[0]-1
            alignment[0].append(seq1[row_idx])
            alignment[1].append('')
            alignment[2].append('-')
        else:
            start[1]=start[1]-1
            alignment[0].append('-')
            alignment[1].append('')
            alignment[2].append(seq2[col_idx])
    align1 = list(reversed(alignment[0]))
    matches = list(reversed(alignment[1]))
    align2 = list(reversed(alignment[2]))
    return (align1, matches, align2)


def calculate_scores(matches, matrix, start):
    """calculate the alignment score and identity score

    matches:
    matrix: list of list; a matrix consisting of a list of list containing
        the best possible alignment score per position.
    start: list; The final end position this is where the traceback needs
        to start [row index, colom index].

    :return: None
    """
    alg_score =matrix[start[0]][start[1]]
    ident_score= matches.count('|')/len(matches)*100
    return (alg_score,ident_score)


def main(seq1,seq2, end_gap_pen, gap_pen):
    """this is the main finction of the script"""
    blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""
    order, blosum = blosum_parser(blosum)
    matrix = initial_matrix(seq1,seq2,end_gap_pen)
    matrix,start,traceback=fil_matrix(seq1,seq2,matrix,gap_pen,order,blosum)
    traceback=completing_traceback(traceback, seq1, seq2)
    align1,matches,align2=traceback_matrix(traceback, start, seq1, seq2)
    alg_score,ident_score= calculate_scores(matches,matrix,start)
    return (align1,align2,alg_score,traceback, ident_score)


if __name__ == "__main__":
    """Here the exercisess are run"""
    print("Question 1")
    #figure 5.9
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    end_gap_pen=-8
    gap_pen=-8
    # NM: I would call the function not main(), but get_alignment()
    align1,align2,alg_score, traceback, ident_score= main(seq1,seq2,
                                                          end_gap_pen, gap_pen)
    print("alignment penalty -8 (figure 5.9):")
    print(*align1, sep=' ')
    print(*align2, sep=" ")
    #figure 5.11
    end_gap_pen = -4
    gap_pen = -4
    align1,align2,alg_score,traceback,ident_score = main(seq1,seq2
                                                         ,end_gap_pen,gap_pen)
    print("alignment penalty -4 (figure 5.11):")
    print(*align1, sep=' ')
    print(*align2, sep=" ")
    #figure 5.12
    end_gap_pen = 0
    gap_pen = -8
    align1,align2,alg_score,traceback, ident_score = main(seq1,seq2,
                                                          end_gap_pen,gap_pen)
    print("alignment end penalty 0, penalty -8 (figure 5.12):")
    print(*align1, sep=' ')
    print(*align2, sep=" ")
    seq3 = ("MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQ"
    "TGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKD"
    "IAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVG"
    "ENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPC"
    "FEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRV"
    "FKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA") #rows
    seq4 = ("MGLLCSRSRHHTEDTDENAQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQASS"
    "DKRKIIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDPAKYTLSSEN"
    "MAIGEKLSEIGARLDYPRLTKDLAEGIETLWNDPAIQETCSRGNELQVPDCTKYLMENLK"
    "RLSDVNYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHL"
    "FEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSIMLFLNKFDIF"
    "EKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYR"
    "TTALDQKLVKKTFKLVDETLRRRNLLEAGLL") #columns
    end_gap_pen = -1
    gap_pen = -5
    print("Question 2:\n", traceback)
    print("Question 3:")
    align1,align2,alg_score,traceback, ident_score =main(seq3,seq4,
                                                         end_gap_pen,gap_pen)
    print("alignment score end penalty 1, penalty 5:",alg_score)
    print("identity score end penalty 1, penalty 5:", ident_score)
    end_gap_pen = -5
    gap_pen = -10
    align1,align2,alg_score,traceback, ident_score =main(seq3,seq4,
                                                         end_gap_pen,gap_pen)
    print("alignment score end penalty 5, penalty 10:", alg_score)
    print("identity score end penalty 5, penalty 10:", ident_score)