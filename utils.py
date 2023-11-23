import os
import multiprocessing as mp
from pydivsufsort import divsufsort, sa_search
from typing import List

def seed_and_extend(
    seq: str,
    bait: str,
    m: int,
    sa = None,
    seed_len: int = 10,
    seqlens = None,
):
    '''
    Applies the seed-and-extend heuristic to search for imperfect matches of the bait
    in the sequence. Returns a list of indices of the matches.

    Arguments:
        seq (str): Sequence to search for the baits in. 
        bait (str): Bait to search for.
        m (int): Allowed Hamming distance for each match.
        sa (suffix array): Suffix array of the string. If not provided, will be constructed.
        seed_len (int): Seed length to use for the heuristic. Every seed-length substring of the 
            bait will be used as seeds.
        seqlens (list[int]): If `seq` was obtained by concatenating sequences, this list 
            should have the lengths of each individual sequence before concatenation. This 
            is because baits should not be able to align around the concatenation points, so
            these values will be used to eliminate those baits.
    '''
    if sa is None:
        # checkpoint = time.time()
        sa = divsufsort(seq)
        # print('Divsufsort took {} seconds.'.format(time.time() - checkpoint))
    if seqlens is None:
        seqlens = [len(seq)]

    lb = len(bait)
    if seed_len > lb:
        raise Exception('Specified seed length is larger than the actual bait.')
    
    final_matches = set()
    checked_indices = set()

    # lookup_time = 0
    # discard_time = 0
    # alignment_time = 0
    # d_picked = 0
    # d_checked = 0
    # d_concat = 0
    # d_negative = 0
    
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
        sum_seqlens = 0
        for seqlen in seqlens:
            discard_matches = set()
            sum_seqlens += seqlen
            for match in seed_matches:
                # find actual start of this alignment and check if it covers a concatenation
                actual_start = match - i
                if actual_start in final_matches: # if we already checked the coverage
                    discard_matches.add(match)     
                    # d_picked += 1
                elif actual_start in checked_indices:
                    discard_matches.add(match)
                    # d_checked += 1
                elif actual_start < 0:
                    discard_matches.add(match)
                    # d_negative += 1
                elif actual_start < sum_seqlens and actual_start > sum_seqlens - lb:
                    discard_matches.add(match)
                    # d_concat += 1
                checked_indices.add(actual_start)
            seed_matches = seed_matches.difference(discard_matches) # VERY INEFFICIENT - change
        # discard_time += time.time() - checkpoint
        # checkpoint = time.time()
        for match in seed_matches:
            actual_start = match - i
            sub = seq[actual_start: actual_start + lb]
            if len(sub) != len(bait):
                print(actual_start, 'wtf happened here')
            if hamming_decide(bait, sub, m):
                final_matches.add(actual_start)
        # alignment_time += time.time() - checkpoint
    # print('Lookup time:', lookup_time)
    # print('Discard time:', discard_time)
    # print('Alignment time:', alignment_time)
    # print('Discarded because picked:', d_picked)
    # print('Discarded because checked:', d_checked)
    # print('Discarded because negative:', d_negative)
    # print('Discarded because concat:', d_concat)
    # print('Found alignments:', len(final_matches))
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
            