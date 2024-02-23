from typing import List
import random

MAPPING = {
    'A': 0,
    'G': 1,
    'T': 2,
    'C': 3,
    'N': 4
}

def generate_randomized_strings(input_string, k, d):
    alphabet = ['A', 'G', 'T', 'C']
    result = []
    
    for _ in range(k):
        # Randomly select d positions to modify
        positions_to_modify = random.sample(range(len(input_string)), d)
        
        # Randomize characters at selected positions
        modified_string = list(input_string)
        for pos in positions_to_modify:
            modified_string[pos] = random.choice(alphabet)
        
        result.append(''.join(modified_string))
    
    return result

def get_max_indices(lst):
    max_value = max(lst)
    max_indices = [index for index, elem in enumerate(lst) if elem == max_value]
    return max_indices

def get_max_from_entries(entries):
    freqs = [entry[2] for entry in entries]
    max_value = max(freqs)
    return [entry for entry in entries if entry[2] == max_value]
    
def wfc_csp(
    strings: List[str],
    alphabet: List[str],
    d: int,
    max_iter: int = 100,
    verbose: bool = False
):
    '''
    Apply the WFC-CSP heuristic to find a Hamming-center for the given set of strings.
    '''
    # grab string length
    length = len(strings[0])
    # initialize empty solution strings ('X' is assumed not to be in the alphabet)
    solution = ['X'] * length
    char_freq = [[0 for _ in range(len(alphabet))] for _ in range(length)]
    for i in range(len(strings)):
        string = strings[i]
        for j in range(length):
            char_freq[j][MAPPING[string[j]]] += 1
    scoreboard = [(l, n, char_freq[l][MAPPING[n]]) for l in range(length) for n in MAPPING.keys()]
    scoreboard = sorted(scoreboard, key = lambda x: x[2], reverse = True)

    distances = [length for _ in range(len(strings))]
    decided = [False] * length
    for _ in range(length):
        farthest = strings[random.choice(get_max_indices(distances))]
        entries = []
        for i in range(len(scoreboard)):
            entry = scoreboard[i]
            if (not decided[entry[0]]) and farthest[entry[0]] == entry[1]:
                entries.append(entry)
        picked_entry = random.choice(get_max_from_entries(entries))
        decided[picked_entry[0]] = True
        solution[picked_entry[0]] = picked_entry[1]
        for i in range(len(strings)):
            if strings[i][picked_entry[0]] == picked_entry[1]:
                distances[i] -= 1
    return ''.join(solution), max(distances)
        
            