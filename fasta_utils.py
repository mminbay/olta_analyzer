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
    repeat_length: Tuple[int, int],
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
        repeat_length (tuple[int, int]): Min and max values for the length of the repeats. The
            repeats will be randomized within this range.
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
        logging.basicConfig(filename = log, level = logging.DEBUG)
        
    BASES = ['A', 'G', 'T', 'C']
    vec = [None] * s_length
    total_populated = 0
    repeats = []
    for i in range(n_repeats):
        l = random.randint(repeat_length[0], repeat_length[1])
        repeats.append(''.join(random.choices(BASES, k = l)))
    if verbose:
        logging.info('Generated the following base repeats:')
        for i, repeat in enumerate(repeats):
            logging.info('{}: {}'.format(i, repeat))
    repeat_coverage = repeat_coverage * s_length
    while total_populated < repeat_coverage:
        i = random.randrange(0, s_length)
        repeat_i = random.randrange(len(repeats))
        repeat = repeats[repeat_i]
        if i + len(repeat) >= s_length:
            continue
        for j in range(len(repeat)):
            if vec[i + j] is not None:
                continue
        if verbose:
            logging.info('Populating index {} with a modification on repeat number {}.'.format(i, repeat_i))
        n_modifs = random.randint(0, repeat_noise)
        modif_locs = random.choices(range(len(repeat)), k = n_modifs)
        for j in range(len(repeat)):
            if j in modif_locs:
                vec[i + j] = random.choice([x for x in BASES if x != repeat[j]])
            else:
                vec[i + j] = repeat[j]
        total_populated = len(vec) - vec.count(None)
        generated_repeat = ''.join(vec[i:i + len(repeat)])
        if verbose:
            logging.info('Randomized {} indices:\n{}'.format(n_modifs, modif_locs))
            logging.info('Populated index with the following:\n{}'.format(generated_repeat))
            logging.info('Current total repeat length: {}'.format(total_populated))
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
    print("Number of sequences read: ", read)
    print("Total length: ", length)
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