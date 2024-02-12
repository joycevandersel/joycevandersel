#!/usr/bin/env python3

"""
Author: Joyce van der Sel

Description: this is a script to create a Hidden Markov Model from a sequence
and to generate sequences from this models

Usage: python [script.py]
"""
# Import statements
import random
import time

# Function definitions

# Background amino acid probabilities
pa = {'A': 0.074, 'C': 0.025, 'D': 0.054, 'E': 0.054, 'F': 0.047, 'G': 0.074,
      'H': 0.026, 'I': 0.068, 'L': 0.099, 'K': 0.058, 'M': 0.025, 'N': 0.045,
      'P': 0.039, 'Q': 0.034, 'R': 0.052, 'S': 0.057, 'T': 0.051, 'V': 0.073,
      'W': 0.013, 'Y': 0.034}


def sample(events):
    """Return a key from dict based on the probabilities

    :param events: dictionary; {key: probability}, probabilities can also be
                    weights.
    :return: str; a random key based on the probabilities.
    """
    pick = random.choices(list(events.keys()), list(events.values()))[0]
    return pick


def parse_fasta(fasta):
    """Parsing a fasta file.

    :param fasta: list; every item of the list is a DNA sequence.
    :return: dictionary; the key of the dictionary is sequence header;
                the value is the sequence.
    """
    sequences = {}
    for seq in fasta:
        seq = list(seq.split('\n'))
        sequences[seq[0]] = ''.join(seq[1:])
    return sequences


def guide_list(seqs):
    """ Creates a list that keeps track of a position is a match state or not.

    :param seqs: list; protein sequence alignments.
    :return: list; containing a 'M' or 'I', depending on if it is a match state
                or not.

    A position is a match state when more than half of the positions in
    the sequences contain a protein.
    """
    pos = list(zip(*seqs))
    guide = []
    for item in pos:
        if item.count('-') > len(item) / 2:
            guide.append('I')
        else:
            guide.append('M')
    return guide


def write_reduced_fasta(fasta, name):
    """Makes output file with reduced alignments.
    
    :param fasta: the key of the dictionary is sequence header;
                    the value is the alignment.
    :return: None
    
    The reduced alignment includes only the match states. A position is a match
    state when more than half of the positions in the sequences contain a 
    protein.
    """
    pos = list(zip(*fasta.values()))
    new_pos = []
    for item in pos:
        if item.count('-') < len(item) / 2:
            new_pos.append(item)
    new_seqs = [''.join(x) for x in list(zip(*new_pos))]
    with open(name, 'w') as file:
        for i in range(len(new_seqs)):
            file.write('>' + list(fasta.keys())[i] + '\n')
            file.write(new_seqs[i] + '\n')


def transition_dictionary(seqs, guide):
    """ Creates a transition dictionary from a match state to the next.

    :param seqs: list; protein sequence alignments.
    :param guide: list; containing a 'M' or 'I', depending on if it is a match
            state or not.
    :return: dictionary; keys are the transitions types. The value is a list
                with a length of the number of match states + 1. Each position
                in the list is an integer with the number of times the type of
                transitions happened from the match state until the next.

    This function tracks the transitions sequence by sequence. It starts in a
    match states at position 0.
    """
    names = ['MM', 'MI', 'MD', 'IM', 'II', 'ID', 'DM', 'DI', 'DD']
    match_states = [0 for pos in range(guide.count('M') + 1)]
    trans_dic = {}
    for name in names:
        trans_dic[name] = match_states.copy()
    for item in seqs:
        state = 'M'
        pos = 0
        for i in range(len(guide)):
            if guide[i] == 'M':
                if item[i] != "-":
                    trans_dic[state + 'M'][pos] += 1
                    state = 'M'
                    pos += 1
                else:
                    trans_dic[state + 'D'][pos] += 1
                    state = 'D'
                    pos += 1
            else:
                if item[i] != '-':
                    trans_dic[state + 'I'][pos] += 1
                    state = 'I'
        trans_dic[state + 'M'][-1] += 1
    return trans_dic


def normalize_transition_dic(trans):
    """Normalizes the occurrence score of transitions

    :param trans: dictionary; keys are the transitions types. The value is a
                list with a length of the number of match states + 1. Each
                position in the list is an integer with the number of times the
                type of transitions happened from the match state until the
                next.
    :return: dictionary; the keys are the transition types. The values is a
                list of probabilities per position in a sequence.
    """
    states = [[], [], []]
    for key, value in trans.items():
        if key.startswith('M'):
            states[0].append(value)
        elif key.startswith('I'):
            states[1].append(value)
        else:
            states[2].append(value)
    total_per_state = []
    for state in states:
        per_pos = list(zip(*state))
        sum_per_pos = []
        for pos in per_pos:
            sum_per_pos.append(sum(pos))
        total_per_state.append(sum_per_pos)
    nor_trans = {}
    for key in trans.keys():
        nor_trans[key] = []
    for key, value in trans.items():
        for i in range(len(value)):
            if key.startswith('M'):
                nor_trans[key].append(0 if value[i] == 0 else round(
                    value[i] / total_per_state[0][i], 2))
            elif key.startswith('I'):
                nor_trans[key].append(0 if value[i] == 0 else round(
                    value[i] / total_per_state[1][i], 2))
            else:
                nor_trans[key].append(0 if value[i] == 0 else round(
                    value[i] / total_per_state[2][i], 2))
    return nor_trans


def gather_match_states(seqs):
    """ Calculates the protein options for every match state.

    :param seqs: list of protein sequence alignments.
    :return: dictionary; the key is match state position in a sequence; the
                value is a list of lists with a sublist containing the
                occurrence of the protein.

    A position is a match state when more than half of the positions in
    the sequences contain a protein.
    """
    match = {}
    m_state = 0
    for x in range(len(max(seqs, key=len))):
        options = {}
        for item in seqs:
            if item[x] in options:
                options[item[x]] += 1
            else:
                options[item[x]] = 1
        if '-' in options and options['-'] < len(
                seqs) / 2 or '-' not in options:
            match[m_state] = [[k, v] for k, v in options.items() if k != '-']
            m_state += 1
    return match


def normalize_match_states(matches):
    """Normalize the score of appearance in sequence.

    :param matches:dictionary; the key is match state position in a sequence;
                    the value is a list of tuples with a tuple containing the
                    occurrence of the protein.
    :return:dictionary; the key is match state position in a sequence; the
                value is a list of lists with a sublist containing the
                probability of occurrence of the protein.

    The probability is rounded at two decimals.
    """
    norm_matches = {}
    for key in matches.keys():
        norm_matches[key] = []
    for k, v in matches.items():
        zipped_list = list(zip(*v))
        total = sum(zipped_list[1])
        for i in range(len(v)):
            score = round(v[i][1] / total, 2)
            norm_matches[k].append([v[i][0], score])
    return norm_matches


def sample_HMM(norm_trans, norm_matches):
    """ Creates a sample based on the HMM

    :param norm_trans: dictionary; the keys are the transition types. The
                    values is a list of probabilities per position in a
                    sequence
    :param matches: dictionary; the key is match state position in a sequence;
                    the value is a list of lists with a sublist containing the
                    probability of occurrence of the protein.
    :return: str; a sequence generated based on the probabilities in the given
                    dictionaries

    This function calls on the function sample, which based on the probability
    chooses a random key to report back
    """
    state = 'M'
    pos = 0
    sample_seq = ""
    while pos != len(list(norm_trans.values())[0]) - 1:
        event = {}
        if state == "M":
            for k, v in norm_trans.items():
                if k.startswith('M'):
                    event[k] = v[pos]
            action = sample(event)
            aa = sample(dict(norm_matches[pos]))
            sample_seq += aa
            state = action[1]
            if state != 'I':
                pos += 1
        elif state == "I":
            for k, v in norm_trans.items():
                if k.startswith('I'):
                    event[k] = v[pos]
            action = sample(event)
            aa = sample(pa)
            sample_seq += aa
            state = action[1]
            if state != 'I':
                pos += 1
        else:
            for k, v in norm_trans.items():
                if k.startswith('D'):
                    event[k] = v[pos]
            action = sample(event)
            sample_seq += '-'
            state = action[1]
            if state != 'I':
                pos += 1
    return sample_seq


def main(infile, name):
    """This is the main function of the script

    :param infile: list; every item is a sequence combined with a header
    :return: tuple of two dictionaries; dictionary one: the key is match state
                position in a sequence; the value is a list of lists with a
                sublist containing the occurrence of the protein. The second
                dictionary is
    """
    fasta = parse_fasta(infile[1:])  # first is empty
    guide = guide_list(list(fasta.values()))
    write_reduced_fasta(fasta, name)
    trans_dic = transition_dictionary(list(fasta.values()), guide)
    norm_trans = normalize_transition_dic(trans_dic)
    protein_bases = gather_match_states(list(fasta.values()))
    norm_matches = normalize_match_states(protein_bases)
    sample_seq = sample_HMM(norm_trans, norm_matches)
    return norm_trans, norm_matches


def print_question_answers(norm_trans, name, norm_matches, name2,
                           norm_trans_large, norm_matches_large):
    print("Question 1: {} match states".format(len(norm_trans) - 1))
    # minus one for the end state
    print("Question 2: the file is called '{}'".format(name))
    print("Question 3:")
    print("transition:", norm_trans)
    print("emission:", norm_matches)
    print('pa:', pa)
    print("Question 4:")
    for x in range(1, 11):
        sample_seq = sample_HMM(norm_trans, norm_matches)
        print('\t', x, sample_seq)
    print("Question 5: \n \t the file is called '{}'".format(name2))
    for x in range(1, 11):
        sample_seq = sample_HMM(norm_trans_large, norm_matches_large)
        print('\t', x, sample_seq)
    # print(norm_trans_large)
    # print(norm_matches_large)


if __name__ == "__main__":
    start_time = time.time()
    infile = open('test.fasta').read().split('>')
    name = 'reduced_test.fasta'
    norm_trans, norm_matches = main(infile, name)
    infile = open('test_large.fasta').read().split('>')
    name2 = 'reduced_large.fasta'
    norm_trans_large, norm_matches_large = main(infile, name2)
    print_question_answers(norm_trans, name, norm_matches, name2,
                           norm_trans_large, norm_matches_large)
    end_time = time.time()
    # print(end_time - start_time)
