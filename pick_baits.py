from cluster_tools import ClusterHelper
from fasta_utils import sequences_to_fasta, sequences_from_fasta
import logging
import multiprocessing as mp
import os
from pydivsufsort import divsufsort
import sys
import time
from tqdm import tqdm
from typing import Set, List, Union, Tuple
from utils import *

'''
    This file is meant to take care of picking baits for a set of sequences.
    The clustering should already be taken care of.
'''
def pick_baits_syotti(
    s: Union[str, List[str]], 
    l: int,
    m: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
) -> Tuple[List[int], List[str]]:
    '''
    Algorithm to pick baits imitating the Syotti method. Unlike the original
    Syotti implementation, reverse complements are optional.

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length.
        m (int): Mismatch allowance.
        seed_length (int): Seed length to be used for the seed-and-extend heuristic.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
    
    Returns:
        ids (list[int]): Indices of picked baits in the concatenated sequence.
        baits (list[str]): Picked baits.
    '''
    start = time.time()
    if log is not None:
        verbose = True
        logging.basicConfig(filename = log, level = logging.DEBUG)
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = compute_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    length = len(s)
    cov = [False] * length
    if verbose:
        logging.info('Initialized bit array with length {}'.format(length))

    final_baits = []
    final_indices = []

    current_seqend = 0
    seqlens_ctr = -1
    for i in tqdm(range(length), desc = 'Covering sequences', unit = 'indices'):
        # if this index is covered, nothing to do
        if cov[i]:
            continue
        # if this index is greater than the current_seqend, shift sequences
        while i >= current_seqend:
            seqlens_ctr += 1
            current_seqend += seqlens[seqlens_ctr]
        if verbose:
            logging.info('Index {} is not yet covered.'.format(i))
        # if this bait covers a concatenation spot, pick the last bait that doesn't
        if i > current_seqend - l:
            bait_start = current_seqend - l
            if verbose:
                logging.info('Bait at {} would cover a concatenation spot at {}'.format(i, current_seqend))
        else:
            bait_start = i
        if verbose:
            logging.info('Aligning bait at {}'.format(bait_start))
        bait = s[bait_start: bait_start + l]
        matches = list(seed_and_extend(s, bait, m, sa, seed_length, seqlens))
        coverages = calculate_coverage(matches, l)
        if verbose:
            logging.info('Bait covers between:\n {}'.format(coverages))
        for c in coverages:
            for j in range(c[0], c[1]):
                cov[j] = True
        final_baits.append(bait)
        final_indices.append(bait_start)
        # reverse complement
        if rc:
            bait = reverse_complement(bait)
            matches = list(seed_and_extend(s, bait, m, sa, seed_length, seqlens))
            coverages = calculate_coverage(matches, l)
            if verbose and len(coverages) > 0:
                logging.info('Reverse complement covers between:\n {}'.format(coverages))
            for c in coverages:
                for j in range(c[0], c[1]):
                    cov[j] = True
    end = time.time()
    print('Remaining uncovered indices: ', cov.count(False))
    print('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    if verbose:
        logging.info('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    return final_indices, final_baits

def pick_baits_syotti_smart(
        s: Union[str, List[str]], 
        l: int,
        m: int,
        seed_length: int = 10,
        rc: bool = False,
        log: str = None,
    ) -> Tuple[List[int], List[str]]:
    '''
    Algorithm to pick baits with a modified syotti method. Instead of committing to
    the first pick that covers an uncovered index, this one will consider all baits
    that cover that index and pick the one with maximum coverage.

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length.
        m (int): Mismatch allowance.
        seed_length (int): Seed length to be used for the seed-and-extend heuristic.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
    
    Returns:
        ids (list[int]): Indices of picked baits in the concatenated sequence.
        baits (list[str]): Picked baits.
    '''
    start = time.time()
    if log is not None:
        verbose = True
        logging.basicConfig(filename = log, level = logging.DEBUG)
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = compute_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    length = len(s)
    cov = [False] * length
    if verbose:
        logging.info('Initialized bit array with length {}'.format(length))

    final_baits = []
    final_indices = []

    current_seqstart = 0
    current_seqend = 0
    seqlens_ctr = -1
    for i in tqdm(range(length), desc = 'Covering sequences', unit = 'indices'):
        # if this index is covered, nothing to do
        if cov[i]:
            continue
        # if this index is greater than the current_seqend, shift sequences
        while i >= current_seqend:
            current_seqstart = current_seqend
            seqlens_ctr += 1
            current_seqend += seqlens[seqlens_ctr]
        if verbose:
            logging.info('Index {} is not yet covered.'.format(i))
        # if this bait covers a concatenation spot, pick the last bait that doesn't
        range_start = i - l + 1
        range_end = i
        if range_start < current_seqstart:
            range_start = current_seqstart
            if verbose:
                logging.info('Baits after {} would cover a concatenation spot at {}'.format(range_start, current_seqstart))
        if range_end > current_seqend - l:
            range_end = current_seqend - l
            if verbose:
                logging.info('Baits before {} would cover a concatenation spot at {}'.format(range_end, current_seqend))
        if verbose:
            logging.info('Considering baits from {} to {}'.format(range_start, range_end))
        max_bait = None
        max_cov = 0
        max_index = -1
        max_covered = None
        for j in range(range_start, range_end + 1):
            bait = s[j: j + l]
            matches = list(seed_and_extend(s, bait, m, sa, seed_length, seqlens))
            if rc:
                matches.extend(list(seed_and_extend(s, reverse_complement(bait), m, sa, seed_length, seqlens)))
            coverages = calculate_coverage(matches, l)
            curr_cov = 0
            for c in coverages:
                curr_cov += cov[c[0]: c[1]].count(False)
            if curr_cov > max_cov:
                max_bait = bait
                max_cov = curr_cov
                max_index = j
                max_covered = coverages
        final_baits.append(max_bait)
        final_indices.append(max_index)
        for c in max_covered:
            for j in range(c[0], c[1]):
                cov[j] = True
    end = time.time()
    print('Remaining uncovered indices: ', cov.count(False))
    print('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    if verbose:
        logging.info('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    return final_indices, final_baits

def pick_baits_vsearch(
        s: Union[str, List[str]], 
        l: int, 
        k: int, 
        dir: str, 
        output: str, 
        log: str = None
    ) -> Set[int]:
    '''
    Algorithm to pick baits. Given a sequence (list) and its VSEARCH clustering
    output of its substrings, outputs a set of baits that cover the entire sequence(s).

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length, or substring length that was used in VSEARCH.
        k (int): Step size that was used in VSEARCH.
        dir (str): Path to the VSEARCH clustering output directory.
        output (str): Output path.
        log (str): If not None, the process will log its progress to the specified path.
    '''
    if log is not None:
        verbose = True
        logging.basicConfig(filename = log, level = logging.DEBUG)
    else:
        verbose = False
    final_baits = set()

    # if multiple sequences were provided, concatenate them in order
    if isinstance(s, list):
        seqlens = [len(seq) for seq in s]
        s = ''.join(s)
    else:
        seqlens = [len(s)]
    
    clusters = ClusterHelper(
        l,
        k,
        seqlens,
        cluster_dir = dir
    )
    
    if verbose:
        logging.info('Created a cluster helper from provided map, with bait length {} and step size {}'.format(l, k))
    
    length = len(s)
    cov = [False] * length
    if verbose:
        logging.info('Initialized bit array with length {}'.format(length))

    prev = None
    indices = {0, length - 1}
    while indices:
        if verbose:
            logging.info('Current endpoints to cover are {}'.format(str(indices)))
        if indices == prev:
            if verbose:
                logging.error('No progress has been made since last interation - killing process')
            raise Exception('No progress has been made - killing process')
        prev = indices
        max_bait = -1
        max_cov = 0
        covered_subs = set()
        checked_baits = set()
        for i in indices: # hepsi yerine bir kismina bakilabilir
            baits = clusters.get_repr(i) # implement
            if verbose:
                logging.info('Index {} is covered by baits {}'.format(str(i), str(baits)))
            for bait in baits:
                # if verbose:
                #     logging.info('Checking bait {}'.format(str(bait)))
                if bait in checked_baits or bait in final_baits:
                    # if verbose:
                    #     logging.info('Already checked or included - skipping bait')
                    continue
                checked_baits.add(bait)
                this_subs = calculate_coverage(clusters.get_subs(bait), l)
                this_cov = 0
                for c in this_subs:
                    this_cov += cov[c[0]: c[1]].count(False)
                # if verbose:
                #     logging.info('Bait {} covers {} previously uncovered cells'.format(str(bait), str(this_cov)))
                if this_cov > max_cov:
                    if verbose:
                        logging.info('Bait {} is the new max coverage bait'.format(str(bait)))
                    max_bait = bait
                    max_cov = this_cov
                    covered_subs = this_subs
        if verbose:
            logging.info('Bait {} is the picked bait'.format(str(max_bait)))
        for c in covered_subs:
            for i in range(c[0], c[1]):
                cov[i] = True
            if verbose:
                logging.info('Covered {} to {}'.format(str(c[0]), str(c[1])))
            indices_remove = set()
            for i in indices:
                if i >= c[0] and i < c[1]:
                    indices_remove.add(i)
            indices = indices - indices_remove
            if c[0] - 1 > 0 and not cov[c[0] - 1]:
                indices.add(c[0] - 1)
            if c[1] < length and not cov[c[1]]:
                indices.add(c[1])
        if verbose:
            logging.info('New indices are {}'.format(str(indices)))
        final_baits.add(max_bait)
    
    final_baits_list = list(final_baits)
    final_baits_str = [s[bait:bait + l] for bait in final_baits_list]
    sequences_to_fasta(final_baits_str, output, final_baits_list)
    return final_baits

def validate_baits(s, l, k, d, baits, dir):
    if isinstance(s, list):
        seqlens = [len(seq) for seq in s]
        s = ''.join(s)
    else:
        seqlens = [len(s)]
    
    clusters = ClusterHelper(
        l,
        k,
        seqlens,
        cluster_dir = dir
    )

    length = len(s)
    cov = [False] * length

    for bait in baits:
        bait_str = s[bait: bait + l]
        subs = clusters.get_subs(bait)
        for sub in subs:
            sub_str = s[sub: sub + l]
            distance = hamming_distance(bait_str, sub_str)
            if distance > d:
                print("WARNING: Bait {} aligns to substring {} with distance {}".format(bait, sub , distance))
                print("Actual bait: {}".format(bait_str))
                print("Aligned substring: {}".format(sub_str))
            else:
                for i in range(sub, sub + l):
                    cov[i] = True
    print("All baits aligned")
    print("Remaining uncovered indices: ", cov.count(False))

if __name__ == '__main__':
    method = sys.argv[1]
    if method == 'syotti':
        spath = sys.argv[2]
        l = int(sys.argv[3])
        m = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        _, sequences = sequences_from_fasta(spath)
        ids, baits = pick_baits_syotti(
            sequences,
            l,
            m,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        sequences_to_fasta(baits, output, ids = ids)
    elif method == 'syotti_s':
        spath = sys.argv[2]
        l = int(sys.argv[3])
        m = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        _, sequences = sequences_from_fasta(spath)
        ids, baits = pick_baits_syotti_smart(
            sequences,
            l,
            m,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        sequences_to_fasta(baits, output, ids = ids)
    else:
        raise Exception('what')
        



    # print('Starting to pick baits...')
    # start = time.perf_counter()

    # res = pick_baits_vsearch(
    #     sequences, 
    #     l, 
    #     k, 
    #     cluster_dir,
    #     output,
    #     log = output[:output.rfind('.')] + '.log'
    # )

    # end = time.perf_counter()
    # print('Picked baits in {} seconds'.format(str(end - start)))
    # print(sorted(list(res)))
    # print(len(res))
    # validate_baits(sequences, l, k, 40, res, cluster_dir)
    # print('Validated baits in {} seconds'.format(str(end - time.perf_counter())))