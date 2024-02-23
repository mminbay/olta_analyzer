import csv
import matplotlib.pyplot as plt
from fasta_utils import sequences_from_fasta
import os
from pydivsufsort import divsufsort
import statistics
import sys
from typing import List, Union, Tuple
from utils import seed_and_extend, calculate_seqlens, calculate_coverage

def calculate_redundancy(
    ids: List[str],
    baits: List[str],
    s: List[str],
    m: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None
):
    '''
    Given a list of sequences and a solution bait set, returns an integer vector that indicates
    how many baits cover an index in the concatenated sequence list.
    Accuracy of results repend on seed length. For consistency, use the same seed length
    provided to the picking algorithm.

    Arguments:
        ids (list[str]): IDs of solution baits.
        baits (list[str]): Solution baits.
        seqs (list[str]): Set of sequences that the solution baits cover.
        m (int): Allowed mismatches for bait binding.
        seed_length (int): Seed length to be used in the seed and extend algorithm.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
    '''
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Calculating coverarage redundancy with provided arguments:\n')
        f.write('m (mismatch allowance) = {}\n'.format(m))
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
    cov = [0] * length
    if verbose:
        f.write('Initialized integer array with length {}\n'.format(length))

    for i, bait in enumerate(baits):
        if verbose:
            f.write('Aligning bait {}.\n'.format(ids[i]))
        l = len(bait)
        matches = list(seed_and_extend(s, bait, m, sa, seed_length, seqlens))
        coverages = calculate_coverage(matches, l)
        if verbose:
            f.write('Bait covers between:\n {}\n'.format(coverages))
        for c in coverages:
            for j in range(c[0], c[1]):
                cov[j] += 1
    if verbose:
        f.write('--------\n')
        f.write('Remaining uncovered indices: {}.\n'.format(cov.count(0)))
        f.write('Total coverage: {}.\n'.format(sum(cov)))
        f.write('Total redundant coverage: {}.\n'.format(sum(cov) - len(cov) + cov.count(0)))
        f.write('Redundancy: {}.\n'.format((sum(cov)/(len(cov) - cov.count(0))) - 1))
        f.write('Mean coverage: {}.\n'.format(statistics.mean(cov)))
        f.write('Standard deviation: {}.\n'.format(statistics.stdev(cov)))
        f.write('Range: {}-{}.\n'.format(min(cov), max(cov)))
        f.close()
    return cov

def calculate_fairness(
    ids: List[str],
    baits: List[str],
    s: List[str],
    m: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None
):
    '''
    Given a list of sequences and a solution bait set, returns a dictionary that indicates
    the coverage of each provided bait.
    Accuracy of results repend on seed length. For consistency, use the same seed length
    provided to the picking algorithm.

    Arguments:
        ids (list[str]): IDs of solution baits.
        baits (list[str]): Solution baits.
        seqs (list[str]): Set of sequences that the solution baits cover.
        m (int): Allowed mismatches for bait binding.
        seed_length (int): Seed length to be used in the seed and extend algorithm.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
    '''
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Calculating fairness with provided arguments:\n')
        f.write('m (mismatch allowance) = {}\n'.format(m))
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

    result = {}
    for id in ids:
        result[id] = 0

    for i, bait in enumerate(baits):
        if verbose:
            f.write('Aligning bait {}.\n'.format(ids[i]))
        l = len(bait)
        matches = list(seed_and_extend(s, bait, m, sa, seed_length, seqlens))
        coverages = calculate_coverage(matches, l)
        if verbose:
            f.write('Bait covers between:\n {}\n'.format(coverages))
        for c in coverages:
            result[ids[i]] += c[1] - c[0]         
    if verbose:
        f.write('--------\n')
        f.write('Total coverage: {}.\n'.format(sum(result.values())))
        f.write('Mean coverage: {}.\n'.format(statistics.mean(result.values())))
        f.write('Standard deviation: {}.\n'.format(statistics.stdev(result.values())))
        f.write('Range: {}-{}.\n'.format(min(result.values()), max(result.values())))
        f.close()
    return result

if __name__ == '__main__':
    option = sys.argv[1]
    if option == 'redundancy':
        bpath = sys.argv[2]
        spath = sys.argv[3]
        m = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        ids, baits = sequences_from_fasta(bpath)
        _, sequences = sequences_from_fasta(spath)
        vec = calculate_redundancy(
            ids,
            baits,
            sequences,
            m,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        plt.plot(vec, '-')
        plt.xlabel('Indices')
        plt.ylabel('Coverage')
        plt.yticks(range(max(vec) + 1))
        plt.title('Redundancy analysis of {} with {}'.format(
            os.path.basename(spath), 
            os.path.basename(bpath)
        ))
        plt.savefig(output)
        
    elif option == 'fairness':
        bpath = sys.argv[2]
        spath = sys.argv[3]
        m = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        ids, baits = sequences_from_fasta(bpath)
        _, sequences = sequences_from_fasta(spath)
        coverages = calculate_fairness(
            ids,
            baits,
            sequences,
            m,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log'
        )
        
        # Write the dictionary to a CSV file
        with open(output, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Bait ID', 'Coverage'])
            for key, value in coverages.items():
                writer.writerow([key, value])
    else:
        print('moe')
    
        
    