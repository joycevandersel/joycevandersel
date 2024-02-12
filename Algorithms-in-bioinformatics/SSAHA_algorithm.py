#!/usr/bin/env python3
"""
Author: Joyce van der Sel
Student number: 1091565

command line: python [script.py] [fasta file] [query file]
(Both need to be in fasta file format)

Implementation of the SSAHA algorithm. Using indexing to find possible
alignment locations.

This script has been optimized to be able to handle multiple longest hits.
Or Use a length threshold of which hit to report back.
"""
# import statements
import re
import time
from sys import argv


# implement your functions here
def build_hash_table(seqs, k):
    """Here a hash table is build.

    :param seqs:list; a list of sequences.
    :param k: int; k-mer length used for the hash table.
    :return: dict; a hash table containing all occurrences of k-mers
    {k-mer: [(sequence,position)].
    """
    hash_table = {}
    for seq_id, seq in enumerate(seqs):
        for i in range(0, len(seq) - k + 1, k):
            if seq[i:i + k] in hash_table:
                hash_table[seq[i:i + k]].append((seq_id + 1, i + 1))
            else:
                hash_table[seq[i:i + k]] = [(seq_id + 1, i + 1)]
    return hash_table


def create_m_list(hash_table, query, k):
    """ Creates a master list based on a query and the hash_table.

    :param hash_table: dict; a hash table containing all occurrences of k-mers
    {k-mer: [(sequence,position)].
    :param query: str; A DNA sequence string.
    :param k: int; k-mer length used to build the master list (must be the same
                as the k used to build the hash table).
    :return: list of tuples; each tuple containing the index, shift and offset.
    """
    m_list = []
    for i in range(len(query) - k + 1):
        if query[i:i + k] in hash_table:
            hits = hash_table[query[i:i + k]]
            [m_list.extend([(pos[0], pos[1] - i, pos[1])]) for pos in hits]
    m_list.sort()
    return m_list


def find_best_hit(m_list, k):
    """Searches through the master list to connect k-mer hits.

    :param m_list:list of tuples; a tuple containing the index, shift and
                    offset.
    :param k: int; k-mer length used in the hash table.
    :return: list of list; each list is a maximum hits (containing longest
                ungapped match). Each sublist consists of sequence index,
                shift, start and stop offset.

    I chose to implement the function so that it can easily be altered to
    not only return the longest stretch(es) but return all above a threshold
    length. See the commented section at the end of the function.
    """
    matches = {}
    for x in range(len(m_list) - 1):
        if m_list[x][0] == m_list[x + 1][0] and m_list[x][1] == m_list[x + 1][
            1] \
                and m_list[x][2] - m_list[x + 1][2] == -k:
            if (m_list[x][0], m_list[x][1]) in matches:
                matches[(m_list[x][0], m_list[x][1])].append(m_list[x + 1][2])
            else:
                matches[(m_list[x][0], m_list[x][1])] = [m_list[x][2],
                                                         m_list[x + 1][2]]
    # For returning max
    if len(matches) == 0:
        return []
    max_value = max(matches.values(), key=len)
    max_keys = [k for k, v in matches.items() if v == max_value]
    top_hit = []
    for item in max_keys:
        top_hit.append([item[0], item[1]])
        top_hit[-1].append(matches[item][0])
        top_hit[-1].append(matches[item][-1])
    # for returning above threshold (here set at 1)
    # threshold = 1
    # for key, value in matches.items():
    # if len(value) > threshold:
    # top_hit.append([key[0], key[1]])
    # top_hit[-1].append(matches[key][0])
    # top_hit[-1].append(matches[key][-1]+k)
    return top_hit


def prep_alignment(seqs, query, match, k):
    """Here the data from the test alignment is prepared for printing

    :param seqs: list; list of reference DNA sequences
    :param query: str; a DNA sequence string.
    :param match: list of list; each list is a maximum hits (containing longest
                ungapped match). Each sublist consists of sequence index,
                shift, start and stop offset.
    :param k: int; k-mer length used to build the hash table.
    :return: list; an alignment [sequence 1, '|' for matches, sequence 2]
    """
    pos = (match[0][2], match[0][3])
    seq_id = match[0][0]
    shift = match[0][1]
    ref = seqs[seq_id - 1]
    r_seq = list(ref[pos[0] - k * 2:pos[1] + k * 2])
    q_seq = list(query[pos[0] - shift:pos[1] + 2].center(len(r_seq)))
    align = ""
    for x in range(len(r_seq)):
        if q_seq[x] == ' ':
            continue
        elif q_seq[x] == r_seq[x]:
            align += '|'
        else:
            align += '.'
    alignment = [''.join(r_seq), align.center(len(r_seq)), ''.join(q_seq)]
    return alignment


def parse_fasta(ara_genome):
    """Parsing a fasta file

    :param ara_genome: list; every item of the list is a DNA sequence
    :return: tuple of a int and dictionary; containing the total length of the
            sequence, the key of the dictionary is sequence header, the value
            is the sequence
    """
    sequences = {}
    tot_len = 0
    for seq in ara_genome:
        seq = list(seq.split('\n'))
        tot_len += len(''.join(seq[1:]))
        rep_seq = re.sub('[^ACGT]', 'A', ''.join(seq[1:]))
        sequences[seq[0]] = rep_seq
    return tot_len, sequences  # last one is empty


def rev_comp(dna):
    """Creating a complementary DNA strand

    :param dna: str; a DNA sequence
    :return: str; a complementary DNA strand
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev = "".join(complement.get(nuc, nuc) for nuc in reversed(dna))
    return rev


def print_question_answers(hash_table, m_list, seq, tot_len, num_keys, align,
                           query_matches, rev_matches, k):
    """

    :param hash_table:dict; a hash table containing all occurrences of k-mers
    {k-mer: [(sequence,position)].
    :param m_list:list of tuples; each tuple containing the index, shift and
                offset.
    :param seq: int; the number of sequences in the genome (chromosomes)
    :param tot_len: int, the combined length of the sequences
    :param num_keys: int; number of keys in the dictionary representing the
                number of k-mers of length k.
    :param align: list; an alignment [sequence 1, '|' for matches, sequence 2]
    :param query_matches: dic; the key is the query id, the value are the
                lists of forward hits containing (index, shift, start offset,
                and the begin position of the offset stop)
    :param rev_matches: dic; the key is the query id, the value are the
                lists of reversed hits containing (index, shift, start offset,
                and the begin position of the offset stop)
    :param k: int; k-mer length used to build the hash table.
    :return: None
    """
    # Question 1
    print('Question 1:')
    print("\t {} {}".format('K-mer', 'Positions'))
    hash_table = dict(sorted(hash_table.items()))
    for key, value in hash_table.items():
        value = ', '.join('({},{})'.format(*el) for el in value)
        print('\t', key, '\t', value, '\n')
    # Question 2
    print('Question 2:')
    print("\t there are {} hits".format(len(m_list)))
    print("\t the first hit is {} \n".format(m_list[0]))
    print("\t the last hit is {} \n".format(m_list[-1]))

    # Question 3
    print('Question 3:')
    print(align[0])
    print(align[1])
    print(align[2])

    # Question 4
    print('\n Question 4:')
    print("\t Contains {} sequence".format(seq))
    print("\t total sequence length is:", tot_len, '\n')

    # Question 5
    print('Question 5:')
    print("\t it contains {} keys \n".format(num_keys))

    # Question 6
    print('Question 6:')
    for z in query_matches.keys():
        print('\t query {} has {} maximal match(es).'
              .format(k, len(query_matches[z])))
        if len(query_matches[k]) != 0:
            for x in range(len(query_matches[z])):
                print('\t It has a match to chromosome {}.'
                      .format(query_matches[z][x][0]))
                print('\t Starts at position: {}, stops at position: {} \n'
                      .format(query_matches[z][x][2],
                              query_matches[z][x][3] + k))
                # The k is to get to the end of the k-mer this is because

    # Optional
    print('\n Optional Question:')
    for k in query_matches.keys():
        for key in rev_matches.keys():
            if k == key and len(query_matches[k]) != 0:
                print('\t query {} has {} forward match(es).'
                      .format(k, len(query_matches[k])))
                for x in range(len(query_matches[k])):
                    print('\t It has a match to chromosome {}.'
                          .format(query_matches[k][x][0]))
                    print('\t Starts at position: {}, stops at position: {} \n'
                          .format(query_matches[k][x][2],
                                  query_matches[k][x][3] + k))
                    # The k is to get to the end of the k-mer this is because
            if k == key and len(rev_matches[k]) != 0:
                print('\t query {} has {} reversed match(es).'
                      .format(k, len(rev_matches[k])))
                for x in range(len(rev_matches[k])):
                    print('\t It has a match to chromosome {}.'
                          .format(rev_matches[k][x][0]))
                    print('\t Starts at position: {}, stops at position {} \n'
                          .format(rev_matches[k][x][2], rev_matches[k][x]
                    [3] + k))
                    # The k is to get to the end of the k-mer this is because


def main():
    """This is the main function of the script"""
    # For the test data
    query = 'TGCAACAT'

    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1, s2, s3]

    k = 2
    hash_table = build_hash_table(seqs, k)
    m_list = create_m_list(hash_table, query, k)
    match = find_best_hit(m_list, k)
    align = prep_alignment(seqs, query, list(match), k)

    # For the Arabidopsis genome
    k = 13
    ara_genome = open(argv[1]).read().split('>')
    tot_len, genome = parse_fasta(ara_genome[1:])  # first one is empty
    ara_hash = build_hash_table(list(genome.values()), k)
    ara_query = open(argv[2]).read().split('>')
    query_len, ara_query = parse_fasta(ara_query[1:])  # first one is empty
    query_matches = {}
    for k, v in ara_query.items():
        ara_m_list = create_m_list(ara_hash, v, k)
        ara_match = find_best_hit(ara_m_list, k)
        query_matches[k] = ara_match
    # Matches for the reversed query
    rev_matches = {}
    for key, value in ara_query.items():
        rev_seq = rev_comp(value)
        ara_m_list = create_m_list(ara_hash, rev_seq, k)
        ara_match = find_best_hit(ara_m_list, k)
        rev_matches[key] = ara_match
    num_keys = len(ara_hash)
    print_question_answers(hash_table, m_list, len(genome), tot_len, num_keys,
                           align, query_matches, rev_matches, k)


if __name__ == "__main__":
    """Perform these commands when called from the command line"""
    start_time = time.time()
    main()
    end_time = time.time()
    print("time:", end_time - start_time)
