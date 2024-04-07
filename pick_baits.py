from cluster_tools import ClusterHelper
from fasta_utils import sequences_to_fasta, sequences_from_fasta
import multiprocessing as mp
import os
from pydivsufsort import divsufsort
import sys
import time
from tqdm import tqdm
from typing import Set, List, Union, Tuple
from utils import *
from wfc_csp import wfc_csp

'''
    This file is meant to take care of picking baits for a set of sequences.
    The clustering should already be taken care of.
'''

def pick_baits_syotti(
    s: Union[str, List[str]], 
    l: int,
    d: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
) -> Tuple[List[int], List[str]]:
    '''
    NEEDS UPDATE
    Algorithm to pick baits imitating the Syotti method. Unlike the original
    Syotti implementation, reverse complements are optional.

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length.
        d (int): Mismatch allowance.
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
        f = open(log, 'w')
        f.write('Picking baits with provided arguments:\n')
        f.write('l (bait length) = {}\n'.format(l))
        f.write('m (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(rc))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    length = len(s)
    cov = [False] * length
    if verbose:
        f.write('Initialized bit array with length {}\n'.format(length))

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
            f.write('Index {} is not yet covered.\n'.format(i))
        # if this bait covers a concatenation spot, pick the last bait that doesn't
        if i > current_seqend - l:
            bait_start = current_seqend - l
            if verbose:
                f.write('Bait at {} would cover a concatenation spot at {}.\n'.format(i, current_seqend))
        else:
            bait_start = i
        if verbose:
            f.write('Aligning bait at {}.\n'.format(bait_start))
        bait = s[bait_start: bait_start + l]
        matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length),
        coverages = calculate_coverage(matches, l)
        if verbose:
            f.write('Bait covers between:\n {}\n'.format(coverages))
        for c in coverages:
            for j in range(c[0], c[1]):
                cov[j] = True
        final_baits.append(bait)
        final_indices.append(bait_start)
        # reverse complement
        if rc:
            bait = reverse_complement(bait)
            matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length)
            coverages = calculate_coverage(matches, l)
            if verbose and len(coverages) > 0:
                f.write('Reverse complement covers between:\n {}\n'.format(coverages))
            for c in coverages:
                for j in range(c[0], c[1]):
                    cov[j] = True
    end = time.time()
    # print('Remaining uncovered indices: ', cov.count(False))
    # print('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    if verbose:
        f.write('Remaining uncovered indices: {}.\n'.format(cov.count(False)))
        f.write('Picked {} baits in {} seconds.\n'.format(len(final_baits), end - start))
        f.close()
    return final_indices, final_baits

def pick_baits_syotti_smart(
        s: Union[str, List[str]], 
        l: int,
        d: int,
        seed_length: int = 10,
        rc: bool = False,
        log: str = None,
) -> Tuple[List[int], List[str]]:
    '''
    NEEDS UPDATE
    Algorithm to pick baits with a modified syotti method. Instead of committing to
    the first pick that covers an uncovered index, this one will consider all baits
    that cover that index and pick the one with maximum coverage.

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length.
        d (int): Mismatch allowance.
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
        f = open(log, 'w')
        f.write('Picking baits with provided arguments:\n')
        f.write('l (bait length) = {}\n'.format(l))
        f.write('m (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(seed_length))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    length = len(s)
    cov = [False] * length
    if verbose:
        f.write('Initialized bit array with length {}.\n'.format(length))

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
            f.write('Index {} is not yet covered.\n'.format(i))
        # if this bait covers a concatenation spot, pick the last bait that doesn't
        range_start = i - l + 1
        range_end = i
        if range_start < current_seqstart:
            range_start = current_seqstart
            if verbose:
                f.write('Baits before {} would cover a concatenation spot at {}.\n'.format(range_start, current_seqstart))
        if range_end > current_seqend - l:
            range_end = current_seqend - l
            if verbose:
                f.write('Baits after {} would cover a concatenation spot at {}.\n'.format(range_end, current_seqend))
        if verbose:
            f.write('Considering baits from {} to {}.\n'.format(range_start, range_end))
        max_bait = None
        max_cov = 0
        max_index = -1
        max_covered = None
        for j in range(range_start, range_end + 1):
            bait = s[j: j + l]
            matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length)
            if rc:
                matches.extend(seed_and_extend(s, reverse_complement(bait), d, sa = sa, seed_length = seed_length))
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
    # print('Remaining uncovered indices: ', cov.count(False))
    # print('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    if verbose:
        f.write('Remaining uncovered indices: {}.\n'.format(cov.count(False)))
        f.write('Picked {} baits in {} seconds.\n'.format(len(final_baits), end - start))
        f.close()
    return final_indices, final_baits

def pick_baits_syotti_wfc(
        s: Union[str, List[str]], 
        l: int,
        d: int,
        seed_length: int = 10,
        rc: bool = False,
        log: str = None,
) -> Tuple[List[int], List[str]]:
    '''
    NEEDS UPDATE
    Algorithm to pick baits with a modified syotti method, first improvement. In 
    addition to considering all baits that cover a position, this one will use 
    WFC-CSP to compute a Hamming center for the coverage areas and realign that 
    to the sequence.

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length.
        d (int): Mismatch allowance.
        seed_length (int): Seed length to be used for the seed-and-extend heuristic.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
    
    Returns:
        baits (list[str]): Picked baits.
    '''
    start = time.time()
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Picking baits with provided arguments:\n')
        f.write('l (bait length) = {}\n'.format(l))
        f.write('m (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(seed_length))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    length = len(s)
    cov = [False] * length
    ignore = initialize_ignore_vector(seqlens, l)
    if verbose:
        f.write('Initialized bit array and ignore array with length {}.\n'.format(length))

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
            f.write('Index {} is not yet covered.\n'.format(i))
        # if this bait covers a concatenation spot, pick the last bait that doesn't
        range_start = i - l + 1
        range_end = i
        if range_start < current_seqstart:
            range_start = current_seqstart
            if verbose:
                f.write('Baits before {} would cover a concatenation spot at {}.\n'.format(range_start, current_seqstart))
        if range_end > current_seqend - l:
            range_end = current_seqend - l
            if verbose:
                f.write('Baits after {} would cover a concatenation spot at {}.\n'.format(range_end, current_seqend))
        if verbose:
            f.write('Considering baits from {} to {}.\n'.format(range_start, range_end))
        max_bait = None
        max_cov = 0
        max_index = 0
        max_covered = None
        for j in range(range_start, range_end + 1):
            bait = s[j: j + l]
            matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length)
            # IMPLEMENT REVERSE COMPLEMENT
            # if rc:
            #     matches.extend(seed_and_extend(s, reverse_complement(bait), d, ignore, sa = sa, seed_length = seed_length))
            # grab the exact strings from the parts 
            uncovered_match_strings = [s[match: match + l] for match in matches if cov[match: match + l].count(False) > 0]
            # compute center
            center, max_dist = wfc_csp(uncovered_match_strings, ['A', 'G', 'T', 'C', 'N'], d)
            # if center is valid (max_dist <= d), compute its coverage. otherwise, use the initial "exact substring" bait
            if max_dist <= d:
                bait = center
                matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length)
                # IMPLEMENT: REVERSE COMPLEMENT
                # if rc:
                #     matches.extend(seed_and_extend(s, reverse_complement(bait), d, ignore, sa = sa, seed_length = seed_length))
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
    # print('Remaining uncovered indices: ', cov.count(False))
    # print('Picked {} baits in {} seconds.'.format(len(final_baits), end - start))
    if verbose:
        f.write('Remaining uncovered indices: {}.\n'.format(cov.count(False)))
        f.write('Picked {} baits in {} seconds.\n'.format(len(final_baits), end - start))
        f.close()
    return final_indices, final_baits

def pick_baits_wfc_iter(
    s: Union[str, List[str]], 
    l: int,
    d: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
) -> Tuple[List[int], List[str]]:
    '''
    Algorithm to pick baits with a modified syotti method, second improvement. WFC-CSP will be
    iteratively applied to covered regions to find better performing baits until no improvement
    can be made.

    Arguments:
        s (str, or list[str]): Sequence(s) to be covered.
        l (int): Bait length.
        d (int): Mismatch allowance.
        seed_length (int): Seed length to be used for the seed-and-extend heuristic.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
    
    Returns:
        baits (list[str]): Picked baits.
    '''
    start = time.time()
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Picking baits with provided arguments:\n')
        f.write('l (bait length) = {}\n'.format(l))
        f.write('m (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(seed_length))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    length = len(s)
    cov = [False] * length
    ignore = initialize_ignore_vector(seqlens, l)
    if verbose:
        f.write('Initialized bit array and ignore array with length {}.\n'.format(length))
    
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
            f.write('---\nIndex {} is not yet covered.\n'.format(i))
        # if this bait covers a concatenation spot, pick the last bait that doesn't
        range_start = i - l + 1
        range_end = i
        if range_start < current_seqstart:
            range_start = current_seqstart
            if verbose:
                f.write('Baits before {} would cover a concatenation spot at {}.\n'.format(range_start, current_seqstart))
        if range_end > current_seqend - l:
            range_end = current_seqend - l
            if verbose:
                f.write('Baits after {} would cover a concatenation spot at {}.\n'.format(range_end, current_seqend))
        if verbose:
            f.write('Considering baits from {} to {}.\n'.format(range_start, range_end))
        max_bait = None
        max_cov = 0
        max_index = 0
        max_covered = None
        for j in range(range_start, range_end + 1):
            if verbose:
                f.write('-\nConsidering exact substring starting at {} as a bait.\n'.format(j))
            bait = s[j: j + l]
            matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length)
            # indices that are covered by current bait. these are temporarily ignored so later iterations
            # do not realign them
            curr_ignore = []
            # IMPLEMENT: reverse complements
            # if rc:
            #     matches.extend(seed_and_extend(s, reverse_complement(bait), d, ignore, sa = sa, seed_length = seed_length))
            for match in matches:
                curr_ignore.append(match)
                ignore[match] = True
            # find the coverage of the current bait
            coverages = calculate_coverage(matches, l)
            curr_iter_cov = 0 
            for c in coverages:
                curr_iter_cov += cov[c[0]: c[1]].count(False)

            curr_max_cov = curr_iter_cov
            iter = 0
            # to make sure we enter the loop
            improvement = 1
            if verbose:
                f.write('Exact substring aligns to {} regions and covers {} indices.\n'.format(
                    len(matches), 
                    curr_max_cov
                ))
            while improvement > 0:
                iter += 1
                # grab the exact strings from the matches of the last iteration's bait
                if verbose:
                    f.write('Calculating center for aligning regions (iteration {})\n'.format(iter))
                uncovered_match_strings = [s[match: match + l] for match in matches if cov[match: match + l].count(False) > 0]
                # compute center for exact strings
                center, max_dist = wfc_csp(uncovered_match_strings, ['A', 'G', 'T', 'C', 'N'], d)
                # if center is invalid, exit loop. 
                if max_dist > d:
                    if verbose:
                        f.write('Calculated center was invalid: Aborting...')
                    break
                # otherwise, set this center as the new bait and compute its coverage
                bait = center
                new_matches = seed_and_extend(s, bait, d, ignore, sa = sa, seed_length = seed_length)
                # IMPLEMENT: reverse complements
                # if rc:
                #     new_matches.extend(seed_and_extend(
                #         s, 
                #         reverse_complement(bait), 
                #         d,
                #         ignore,
                #         sa = sa,
                #     ))
                for match in new_matches:
                    curr_ignore.append(match)
                    ignore[match] = True
                matches.extend(new_matches)
                coverages = calculate_coverage(matches, l)
                curr_iter_cov = 0
                for c in coverages:
                    curr_iter_cov += cov[c[0]: c[1]].count(False)
                # measure improvement compared to last iteration's bait
                improvement = curr_iter_cov - curr_max_cov
                if verbose:
                    f.write('Iteration {} center aligns to {} regions and covers {} indices.\n'.format(
                        iter,
                        len(matches), 
                        curr_iter_cov
                    ))
                    f.write('{} more indices are covered compared to last iteration.\n'.format(improvement))
                curr_max_cov = curr_iter_cov
            if curr_max_cov > max_cov:
                max_bait = bait
                max_cov = curr_max_cov
                max_index = j
                max_covered = coverages
                if verbose:
                    f.write('Index {}, iteration {} is the new best bait with coverage {}\n'.format(
                        j, iter, max_cov
                    ))
            # revert the indices that were ignored for this iteration
            for index in curr_ignore:
                ignore[index] = False
        final_baits.append(max_bait)
        final_indices.append(max_index)
        for c in max_covered:
            update_ignore_vector(cov, seqlens, ignore, c, l)
            for j in range(c[0], c[1]):
                cov[j] = True
        if verbose:
            f.write('---\nPicked bait covers: {}\n'.format(max_covered))
        if not cov[i]:
            raise Exception('Index {} was somehow not covered\n'.format(i))
    end = time.time()
    if verbose:
        f.write('Remaining uncovered indices: {}.\n'.format(cov.count(False)))
        f.write('Picked {} baits in {} seconds.\n'.format(len(final_baits), end - start))
        f.close()
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
        f = open(log, 'w')
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
        f.write('Created a cluster helper from provided map, with bait length {} and step size {}.\n'.format(l, k))
    
    length = len(s)
    cov = [False] * length
    if verbose:
        f.write('Initialized bit array with length {}.\n'.format(length))

    prev = None
    indices = {0, length - 1}
    while indices:
        if verbose:
            f.write('Current endpoints to cover are {}.\n'.format(str(indices)))
        if indices == prev:
            if verbose:
                f.write('No progress has been made since last interation - killing process.\n')
                f.close()
            raise Exception('No progress has been made - killing process')
        prev = indices
        max_bait = -1
        max_cov = 0
        covered_subs = set()
        checked_baits = set()
        for i in indices: # hepsi yerine bir kismina bakilabilir
            baits = clusters.get_repr(i) # implement
            if verbose:
                f.write('Index {} is covered by baits {}'.format(str(i), str(baits)))
            for bait in baits:
                # if verbose:
                #     f.write('Checking bait {}.\n'.format(str(bait)))
                if bait in checked_baits or bait in final_baits:
                    # if verbose:
                    #     f.write('Already checked or included - skipping bait.\n')
                    continue
                checked_baits.add(bait)
                this_subs = calculate_coverage(clusters.get_subs(bait), l)
                this_cov = 0
                for c in this_subs:
                    this_cov += cov[c[0]: c[1]].count(False)
                # if verbose:
                #     f.write('Bait {} covers {} previously uncovered cells.\n'.format(str(bait), str(this_cov)))
                if this_cov > max_cov:
                    if verbose:
                        f.write('Bait {} is the new max coverage bait.\n'.format(str(bait)))
                    max_bait = bait
                    max_cov = this_cov
                    covered_subs = this_subs
        if verbose:
            f.write('Bait {} is the picked bait.\n'.format(str(max_bait)))
        for c in covered_subs:
            for i in range(c[0], c[1]):
                cov[i] = True
            if verbose:
                f.write('Covered {} to {}.\n'.format(str(c[0]), str(c[1])))
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
            f.write('New indices are {}.\n'.format(str(indices)))
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
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        _, sequences = sequences_from_fasta(spath)
        ids, baits = pick_baits_syotti(
            sequences,
            l,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        sequences_to_fasta(baits, output, ids = ids)
    elif method == 'syotti_smart':
        spath = sys.argv[2]
        l = int(sys.argv[3])
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        _, sequences = sequences_from_fasta(spath)
        ids, baits = pick_baits_syotti_smart(
            sequences,
            l,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        sequences_to_fasta(baits, output, ids = ids)
    elif method == 'syotti_wfc':
        spath = sys.argv[2]
        l = int(sys.argv[3])
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        _, sequences = sequences_from_fasta(spath)
        ids, baits = pick_baits_syotti_wfc(
            sequences,
            l,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        sequences_to_fasta(baits, output, ids = ids)
    elif method == 'wfc_iter':
        spath = sys.argv[2]
        l = int(sys.argv[3])
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        _, sequences = sequences_from_fasta(spath)
        ids, baits = pick_baits_wfc_iter(
            sequences,
            l,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        print(ids[:5])
        print(baits[:5])
        sequences_to_fasta(baits, output, ids = ids)
    else:
        raise Exception('what')
        