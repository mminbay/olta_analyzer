from Bio import SeqIO
import logging
import math
import os
import random
from typing import List, Tuple, Union
from utils import hamming_distance

def synthesize_sequence(
    s_length: int,
    n_repeats: int,
    repeat_length: Union[int, Tuple[int, int]],
    repeat_coverage: float,
    repeat_noise: int,
    log: str = None
) -> List[str]:
    '''
    Generate a randomized nucleotide sequence with the given parameters.

    Arguments:
        s_length (int): Length of the sequence to be generated.
        n_repeats (int): Number of unique repeats. These repeats will be randomly
            scattered around the generated sequence.
        repeat_length (int or tuple[int, int]): Value, or min and max values for the length of the repeats.
        If a tuple if provided, the repeats will be randomized within this range.
        repeat_coverage ([0, 1]): The proportion of the finals sequence that should be populated
            with repeats. Note that this is a lower bound: The randomized regions may have repeats
            by chance.
        repeat_noise (int): Maximum Hamming distance to randomize repeats with. Once a repeat is
            decided, the sequence will be populated with its variants that are at most this distance
            away from that repeat.
        log (str): If provided, will output a log of the process to the provided path.
    '''
    verbose = False
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Generating baits with following arguments:\n')
        f.write('L = {}\n'.format(s_length))
        f.write('RN (repeat number) = {}\n'.format(n_repeats))
        f.write('RL (repeat length) = {}\n'.format(repeat_length))
        f.write('RC (repeat coverage) = {}\n'.format(repeat_coverage))
        f.write('RE (repeat error) = {}\n'.format(repeat_noise))
        
    BASES = ['A', 'G', 'T', 'C']
    vec = [None] * s_length
    total_populated = 0
    repeats_pools = {}
    if verbose:
        repeats_locations = {}
    for i in range(n_repeats):
        l = repeat_length
        if type(repeat_length) is tuple:
            l = random.randint(repeat_length[0], repeat_length[1])
        repeat = ''.join(random.choices(BASES, k = l))
        while repeat in repeats_pools.keys():
            repeat = ''.join(random.choices(BASES, k = l))
        repeats_pools[repeat] = set(range(s_length - l + 1))
        if verbose:
            repeats_locations[repeat] = []
    repeat_coverage = repeat_coverage * s_length
    repeats = list(repeats_pools.keys())
    repeats_ctr = 0
    while total_populated < repeat_coverage:
        repeat = repeats[repeats_ctr]
        i = -1
        if repeat_coverage == s_length:
            i = total_populated # if entire sequence is repeats, fill the next index
        elif len(repeats_pools[repeat]) == 0:
            if verbose:
                f.write('WARNING: {} has no more indices it can be planted in.\n'.format(repeat))
            repeats.remove(repeat)
            n_repeats -= 1
            repeats_ctr = repeats_ctr % n_repeats
            continue
        else:
            i = random.choice(list(repeats_pools[repeat]))
        if i + len(repeat) > s_length:
            if repeat_coverage == s_length:
                vec.extend([None] * (i + len(repeat) - s_length)) # extend vector to accomodate this plant
                if verbose:
                    f.write('Extended sequence to {} base pairs to add final repeat.'.format(len(vec)))
            else:
                continue
        for j in range(len(repeat)):
            if vec[i + j] is not None:
                raise Exception('this shouldnt happen')
        if verbose:
            repeats_locations[repeat].append(i)
        n_modifs = random.randint(0, repeat_noise)
        modif_locs = random.choices(range(len(repeat)), k = n_modifs)
        for j in range(len(repeat)):
            if j in modif_locs:
                vec[i + j] = random.choice([x for x in BASES if x != repeat[j]])
            else:
                vec[i + j] = repeat[j]
        for other in repeats:
            for j in range(i - len(other) + 1, i + len(repeat)):
                if j in repeats_pools[other]:
                    repeats_pools[other].discard(j)
        repeats_ctr = (repeats_ctr + 1) % n_repeats
        total_populated = len(vec) - vec.count(None)
    if verbose:
        f.write('Total covered base pairs: {}\n'.format(total_populated))
        for k, v in repeats_locations.items():
            if len(v) == 0 or len(v) == 1:
                print('WARNING: Some repeats were planted < 2 times')
                f.write('WARNING: Some repeats were planted < 2 times\n')
                break
        f.write('Repeats and the locations they were planted in:\n')
        for k, v in repeats_locations.items():
            f.write(str(k) +  ': ' + str(v) + '\n')
    for i in range(s_length):
        if vec[i] is None:
            vec[i] = random.choice(BASES)
    return ''.join(vec)
        
def sequences_from_fasta(
    file: str,
    n: int = None,
    min_length: int = None,
    indices: List[int] = None,
    start: int = 0,
    slice_seqs: bool = False,
) -> Tuple[List[str], List[str]]:
    '''
    Read a fasta file and return a list of sequences based on the provided arguments.

    Arguments:
        file (str): Path to fasta file
        n (int): Number of sequences to be read. If not None, reading will stop after
            this many sequences have been read.
        min_length (int): Minimum total length of sequences to be read. If not None, 
            reading will stop after total length exceeds this number.
        indices (list(int)): Indices of sequences to be read. Other sequences will be skipped.
        start (int): Index of first sequence to be read. Sequences that appear before this 
            will be skipped.
        slice_seqs (bool): If True, will return a list with a single sequences by concatenating
            all sequences that were read. Will also truncate the sequence to exactly match the
            min_length argument, if it was provided.
    '''
    ids = []
    sequences = []
    i = 0
    read = 0
    length = 0
    for record in SeqIO.parse(file, 'fasta'):
        if i < start:
            i += 1
            continue
        if indices is not None:
            if i not in indices:
                i += 1
                continue
        sequences.append(str(record.seq).upper())
        ids.append(str(record.id))
        length += len(str(record.seq))
        read += 1
        i += 1
        if n is not None:
            if read >= n:
                break
        if min_length is not None:
            if length >= min_length:
                break
    if slice_seqs:
        seq = ''.join(sequences)
        if min_length is not None:
            seq = seq[:min_length]
        sequences = [seq]
        length = len(sequences[0])
    return ids, sequences

def sequences_to_fasta(
    sequences: Union[str, List[str]], 
    output_file: str,
    ids: List[str] = None
):
    '''
    Write a list of sequences to a FASTA file with auto-generated sequence IDs.

    Args:
        sequences (list[str] | str): Sequence, or list of sequences.
        output_file (str): The name of the output FASTA file.
        ids (list[str]): IDs of the sequences in the outputted fasta file. If None,
            IDs will be sequential.
    '''
    if type(sequences) == str:
        sequences = [sequences]
    with open(output_file, 'w') as fasta_file:
        for i, sequence in enumerate(sequences):
            sequence_id = f'seq{i+1}'  # Generate a unique sequence ID
            if ids is not None: 
                sequence_id = ids[i]
            fasta_file.write(f'>{sequence_id}\n')
            fasta_file.write(f'{sequence}\n')