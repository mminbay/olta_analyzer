import sys
import os
import math
from Bio import SeqIO
import logging
from fasta_utils import sequences_to_fasta, sequences_from_fasta
from typing import List

'''
This file takes care of extracting sequences from a larger .fasta file,
and making substrings from the 
'''
def substrings_from_sequences(
    sequences: List[str],
    l: int = 120,
    k: int = 1
):
    '''
    From a list of sequences, return a list of every l-length sequence 
    with k-apart starting indices. The last l-length substring of every 
    sequence is always included, so the entire sequence is covered.

    Arguments:
        sequences (list[str]): Sequences to generate substrings from.
        l (int): Length of substrings to generate.
        k (int): Distance between starting indices of substrings to 
            generate.
    '''
    substrings = []
    indices = []
    shift = 0 # when sequences are concatenated, the i'th substring of the second sequences is actually the i + len(s1)'th substring of the concatenated sequence
    for seq in sequences:
        length = len(seq)
        i = 0
        while i < length - l + 1:
            substrings.append(seq[i:i + l])
            indices.append(i + shift)
            i += k
        # notice if k doesn't divide length - l + 1, there might be uncovered base pairs at the end
        if i - k != length - l: # check if previous step didn't cover the last l indices exactly
            substrings.append(seq[length - l: length])
            indices.append(length - l + shift)
        shift += length
    return substrings, indices

if __name__ == '__main__':
    spath = sys.argv[1]
    l = int(sys.argv[2])
    k = int(sys.argv[3])
    output_path = sys.argv[4]
    min_length = int(sys.argv[5])
    _, s = sequences_from_fasta(
        spath,
        n = None,   
        min_length = min_length,
        ids = None,
        start = 0,
        slice_seqs = True
    )
    write_sequences = True
    if write_sequences:
        sequences_to_fasta(s, output_path[:output_path.rfind('.')] + '_sequences.fasta')
    ss, i = substrings_from_sequences(s, l, k)
    sequences_to_fasta(ss, output_path, ids = i)
    
    
    