import os
import multiprocessing as mp
from pydivsufsort import divsufsort, sa_search
from typing import List, Tuple

def update_ignores(
    cov: List[bool],
    seqlens: List[int],
    ignore: List[bool],
    region: Tuple[int, int],
    l: int
):
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

def seed_and_extend(
    seq: str,
    bait: str,
    d: int,
    ignore: List[bool],
    sa = None,
    seed_len: int = 10,
):
    '''
    Applies the seed-and-extend heuristic to search for imperfect matches of the bait
    in the sequence. Returns a list of indices of the matches.

    Assumes that concatenation spots are correctly marked in the ignore vector.

    Arguments:
        seq (str): Sequence to search for the baits in. 
        bait (str): Bait to search for.
        d (int): Allowed Hamming distance for each match.
        sa (suffix array): Suffix array of the string. If not provided, will be constructed.
        seed_len (int): Seed length to use for the heuristic. Every seed-length substring of the 
        ignore (list[boolean]): Indices that have True entries in this list will not be
            aligned. This should be used to mark indices that would cover concatenation spots or 
            already covered s.t. another alignment would not increase the coverage.
            bait will be used as seeds.
        
    '''
    if sa is None:
        sa = divsufsort(seq)

    if ignore is None:
        ignore = [False] * len(seq)

    # the indices of the ignore list that are modified within the scope of the function.
    # these indices should be reverted at the end of the process
    modified_indices = []

    lb = len(bait)
    if seed_len > lb:
        raise Exception('Specified seed length is larger than the actual bait.')
    
    # checked_indices = set()

    # lookup_time = 0
    # discard_time = 0
    # alignment_time = 0
    # d_checked = 0
    # d_concat = 0
    # d_negative = 0
    
    final_matches = []
    # for i in tqdm(range(lb - seed_len + 1), desc = 'Considering seeds', unit = 'seed'):
    for i in range(lb - seed_len + 1):
        # checkpoint = time.time()
        seed = bait[i: i + seed_len]
        seed_matches = sa_search(seq, sa, seed)
        # lookup_time += time.time() - checkpoint
        if seed_matches[1] is None:
            continue
        seed_matches = set(sa[seed_matches[1]: seed_matches[0] + seed_matches[1]])
        # checkpoint = time.time()
        # seq_start = 0
        # seq_end = 0
        # eliminate alignments that cover concatenation points
        # for seqlen in seqlens:
        #     discard_matches = set()
        #     seq_start = seq_end
        #     seq_end += seqlen
        #     for match in seed_matches:
        #         # find actual start of this alignment and check if it covers a concatenation
        #         actual_start = match - i
        #         if ignore[actual_start]:
        #             continue
        #         if actual_start in checked_indices:
        #             discard_matches.add(match)
        #             print('discarded bc already checked')
        #             d_checked += 1
        #         elif actual_start < 0:
        #             continue
        #         elif actual_start < 0:
        #             discard_matches.add(match)
        #             print('discarded bc negative')
        #             d_negative += 1
        #         elif actual_start < seq_end and actual_start > seq_end - lb:
        #             modified_indices.append(actual_start)
        #             ignore[actual_start] = True
        #         elif actual_start < seq_end and actual_start > seq_end - lb:
        #             discard_matches.add(match)
        #             print('discarded bc concat spot')
        #             d_concat += 1
        #     seed_matches = seed_matches.difference(discard_matches) # VERY INEFFICIENT - change
        # discard_time += time.time() - checkpoint
        # checkpoint = time.time()
        for match in seed_matches:
            actual_start = match - i
            if actual_start < 0:
                continue
            if ignore[actual_start]:
                continue
            ignore[actual_start] = True # to make sure we don't check this index again
            modified_indices.append(actual_start) # to make sure we revert it to False at the end
            # checked_indices.add(actual_start)
            sub = seq[actual_start: actual_start + lb]
            if len(sub) != len(bait):
                print(actual_start, 'what happened here')
            if hamming_decide(bait, sub, d):
                final_matches.append(actual_start)
        # alignment_time += time.time() - checkpoint
    # print('Lookup time:', lookup_time)
    # print('Discard time:', discard_time)
    # print('Alignment time:', alignment_time)
    # print('Discarded because checked:', d_checked)
    # print('Discarded because negative:', d_negative)
    # print('Discarded because concat:', d_concat)
    # print('Found alignments:', len(final_matches))
    for i in modified_indices:
        ignore[i] = False
    # if len(final_matches) == 0:
    #     raise Exception('how is this possible')
    for m in final_matches:
        if ignore[m]:
            raise Exception("why return ignor")
    return final_matches

def hamming_distance(str1, str2):
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

def hamming_decide(str1, str2, m):
    '''
    '''
    distance = 0
    for i in range(len(str1)):
        if str1[i] == 'N': # character N does not match with anything
            distance += 1 
        elif str1[i] != str2[i]:
            distance += 1
        if distance > m:
            return False
    return True

def calculate_seqlens(seqs: List[str]):
    return [len(seq) for seq in seqs]

def calculate_coverage(subs, l):
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

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    # Create a list of the complement of each nucleotide in the reversed sequence
    sequence = sequence.upper()
    reverse_complement_list = [complement_dict[base] for base in reversed(sequence)]

    # Join the list into a string
    reverse_complement_str = ''.join(reverse_complement_list)

    return reverse_complement_str
            