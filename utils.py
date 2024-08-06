import ctypes
import multiprocessing as mp
from multiprocessing import Pool
import os
from pydivsufsort import divsufsort, sa_search
from typing import List, Tuple, Set

def initialize_ignore_vector(
    seqlens: List[int],
    l: int,
) -> List[bool]:
    '''
    Given a list of sequence lengths, create an ignore vector where the indices that would
    cover concatenation spots are marked. Concatenation spots cannot be aligned to with baits, 
    so they never need to be checked for alignments and can be ignored.
    '''
    ignore = [False] * sum(seqlens)
    seqlen_sum = 0
    for seqlen in seqlens:
        seqlen_sum += seqlen
        for i in range(seqlen_sum - l + 1, seqlen_sum):
            ignore[i] = True
    return ignore

def update_ignore_vector(
    cov: List[bool],
    seqlens: List[int],
    ignore: List[bool],
    region: Tuple[int, int],
    l: int
) -> None:
    '''
    Given a coverage vector, an ignore vector, and a covered region, find the indices that
    no longer need to be checked for alignments and update the ignore vector accordingly. An
    index can be ignored if the index and the l - 1 indices following it are already covered,
    provided that none of those indices are concatenation spots.
    '''
    total_seqlen = len(cov)
    current_seqstart = 0
    current_seqend = 0
    seqlens_ctr = -1

    check_start = region[0] - l + 1
    if check_start < 0:
        check_start = 0
    check_end = region[1] + l
    if check_end > total_seqlen:
        check_end = total_seqlen

    streak = 0 # how many indices in a row are already covered.
    for i in range(check_start, check_end):
        while i >= current_seqend:
            current_seqstart = current_seqend
            seqlens_ctr += 1
            current_seqend += seqlens[seqlens_ctr]
            streak = 0 # if we cross a concatenation spot, streak is back to zero
        if cov[i]:
            streak += 1
        else:
            streak = 0
        if streak >= l:
            ignore[i - l + 1] = True

def naive_wrapper(
    i: int,
    chars_per_core: int,
    seq: str,
    bait: str,
    d: int,
    ignore: List[bool]
):
    l_b = len(bait)
    matches = []
    for j in range(chars_per_core):
        if ignore[j]:
            continue
        if hamming_decide(bait, seq[j: j + l_b], d):
            curr = j + (i * chars_per_core)
            matches.append(curr)
    return matches
    
def naive_alignment(
    seq: str,
    bait: str,
    d: int,
    ignore: List[bool],
    **kwargs
) -> List[int]:
    '''
    Compute coverages in the naive way.

    Arguments:
        seq (str): Sequence to search for the bait in. 
        bait (str): Bait to search for.
        d (int): Allowed Hamming distance for each match.
        ignore (list[boolean]): Indices that have True entries in this list will not be
            aligned. This should be used to mark indices that would cover concatenation spots or 
            already covered s.t. another alignment would not increase the coverage.
        **kwargs:
            num_cores (int): If greater than 1, will parallelize to provided number of cores.
    '''
    num_cores = 1
    if 'num_cores' in kwargs.keys():
        num_cores = kwargs['num_cores']
        
    l_s = len(seq)
    l_b = len(bait)
    chars_per_core = l_s // num_cores

    individual_args = [(
        i,
        chars_per_core,
        seq[i * chars_per_core: (i + 1) * chars_per_core + l_b - 1],
        bait,
        d,
        ignore[i * chars_per_core: (i + 1) * chars_per_core]
    ) for i in range(num_cores)]

    with Pool(num_cores) as pool:
        results = pool.starmap(naive_wrapper, individual_args)

    return [match for result in results for match in result]

def seed_and_extend(
    seq: str,
    bait: str,
    d: int,
    ignore: List[bool],
    sa = None,
    seed_length: int = 10,
) -> List[int]:
    '''
    Applies the seed-and-extend heuristic to search for imperfect matches of the bait
    in the sequence. Returns a list of indices of the matches.

    Assumes that concatenation spots are correctly marked in the ignore vector.

    Arguments:
        seq (str): Sequence to search for the baits in. 
        bait (str): Bait to align.
        d (int): Allowed Hamming distance for each match.
        ignore (list[boolean]): Indices that have True entries in this list will not be
            aligned. This should be used to mark indices that would cover concatenation spots or 
            already covered s.t. another alignment would not increase the coverage.
        sa (suffix array): Suffix array of the string. If not provided, will be constructed.
        seed_length (int): Seed length to use for the heuristic.
    '''
    if sa is None:
        sa = divsufsort(seq)
        
    modified_indices = []

    lb = len(bait)
    if seed_length > lb:
        raise Exception('Specified seed length is larger than the actual bait.')
    
    final_matches = []
    for i in range(lb - seed_length + 1):
        seed = bait[i: i + seed_length]
        seed_matches = sa_search(seq, sa, seed)
        if seed_matches[1] is None:
            continue
        seed_matches = set(sa[seed_matches[1]: seed_matches[0] + seed_matches[1]])
        for match in seed_matches:
            actual_start = match - i
            if actual_start < 0:
                continue
            if ignore[actual_start]:
                continue
            ignore[actual_start] = True # to make sure we don't check this index again
            modified_indices.append(actual_start) # to make sure we revert it to False at the end
            sub = seq[actual_start: actual_start + lb]
            if len(sub) != len(bait):
                print(actual_start, 'what happened here')
            if hamming_decide(bait, sub, d):
                final_matches.append(actual_start)
    for i in modified_indices: # revert the indices ignored for this iteration
        ignore[i] = False

    return final_matches

def hamming_distance(str1, str2) -> int:
    if len(str1) != len(str2):
        raise ValueError("Input strings must have the same length")

    distance = 0
    str1 = str1.upper()
    str2 = str2.upper()
    for i in range(len(str1)):
        if str1[i] == 'N': # character N does not match with anything
            distance += 1 
        elif str1[i] != str2[i]:
            distance += 1
    return distance

def hamming_decide(str1, str2, d) -> bool:
    '''
    '''
    distance = 0
    for i in range(len(str1)):
        if str1[i] == 'N': # character N does not match with anything
            distance += 1 
        elif str1[i] != str2[i]:
            distance += 1
        if distance > d:
            return False
    return True

def calculate_seqlens(seqs: List[str]) -> List[int]:
    return [len(seq) for seq in seqs]

def calculate_coverage(subs: List[int], l: int) -> List[Tuple[int, int]]:
    '''
    Given a list of substring starting indices and substring length, return an actual list of starting and ending indices for the coverage of these substrings.
    Example: suppose we have subs = [5, 10, 15, 45], l = 10. These cover from [5, 25] and [45, 55]
    '''
    if len(subs) == 0:
        return []
    subs = sorted(subs)
    results = []
    curr_start = subs[0]
    curr_end = subs[0] + l
    for sub in subs[1:]:
        if sub <= curr_end: # if this index is covered by the last started interval, extend the interval
            curr_end = sub + l
        else: # otherwise, end last one and start new one
            results.append((curr_start, curr_end))
            curr_start = sub
            curr_end = sub + l
    results.append((curr_start, curr_end)) # end last interval
    return results

def reverse_complement(sequence) -> str:
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    # Create a list of the complement of each nucleotide in the reversed sequence
    sequence = sequence.upper()
    reverse_complement_list = [complement_dict[base] for base in reversed(sequence)]

    # Join the list into a string
    reverse_complement_str = ''.join(reverse_complement_list)

    return reverse_complement_str
            