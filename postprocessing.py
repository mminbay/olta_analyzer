import csv
import matplotlib.pyplot as plt
from fasta_utils import read_sequences
import numpy as np
import os
from pydivsufsort import divsufsort
import statistics
import sys
from typing import List, Union, Tuple
from tqdm import tqdm
import utils
from multiprocessing import Pool

class Polynomial:
    '''
    Polynomial class for finding coverage probabilities and calculating expected values.
    '''
    def __init__(self, coefficients):
        # coefficients should be a list representing the coefficients of the polynomial
        self.coefficients = coefficients

    def __repr__(self):
        # Returns a string representation of the polynomial
        return "Polynomial({})".format(self.coefficients)

    def __add__(self, other):
        # Addition of two polynomials
        # Pad the coefficients with zeros to make their lengths equal
        padded_self = self.coefficients + [0] * (len(other.coefficients) - len(self.coefficients))
        padded_other = other.coefficients + [0] * (len(self.coefficients) - len(other.coefficients))
        # Add the coefficients element-wise
        result_coefficients = [a + b for a, b in zip(padded_self, padded_other)]
        return Polynomial(result_coefficients)

    def __mul__(self, other):
        # Multiplication of two polynomials
        result_degree = len(self.coefficients) + len(other.coefficients) - 1
        result_coefficients = [0] * result_degree
        for i in range(len(self.coefficients)):
            for j in range(len(other.coefficients)):
                result_coefficients[i + j] += self.coefficients[i] * other.coefficients[j]
        return Polynomial(result_coefficients)

def verify_baits(
    ids: List[str],
    baits: List[str],
    s: Union[str, List[str]],
    d: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
    exact: bool = False,
    num_cores = 1
) -> int:
    '''
    Arguments:
        ids (list[str]): IDs of solution baits.
        baits (list[str]): Solution baits.
        s (list[str]): Set of sequences that the solution baits cover.
        d (int): Allowed mismatches for bait binding.
        seed_length (int): Seed length to be used in the seed and extend algorithm.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
        exact (bool): If True, will use naive alignment to align baits instead of seed-and-extend.    
        num_cores (int): If greater than 1, will parallelize naive alignment.
    '''
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Verifying baits with provided arguments:\n')
        f.write('d (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(rc))
        f.write('Exact matching = {}\n'.format(exact))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = utils.calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]
        
    sa = divsufsort(s)
    alignment_function = utils.naive_alignment # for using exact matching
    if not exact:
        alignment_function = utils.seed_and_extend # otherwise
    
    length = len(s)
    cov = [False] * length
    l = len(baits[0])
    ignore = utils.initialize_ignore_vector(seqlens, l)
    if verbose:
        f.write('Initialized integer array and ignore vector with length {}\n'.format(length))
    for id, bait in zip(ids, baits):
        if verbose:
            f.write('Aligning bait {}.\n'.format(id))
        matches = alignment_function(s, bait, d, ignore, sa = sa, num_cores = num_cores)
        coverages = utils.calculate_coverage(matches, l)
        if verbose:
            f.write('Bait covers between:\n {}\n'.format(coverages))
        for c in coverages:
            for j in range(c[0], c[1]):
                cov[j] = True
    if verbose:
        f.write('--------\n')
        f.write('Remaining uncovered indices: {}.\n'.format(cov.count(False)))
        f.write(str([i for i in range(len(cov)) if not cov[i]]))
        f.close()
    return [i for i in range(len(cov)) if not cov[i]]

def calculate_redundancy(
    ids: List[str],
    baits: List[str],
    s: Union[str, List[str]],
    d: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
    exact: bool = False,
    num_cores = 1
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
        d (int): Allowed mismatches for bait binding.
        seed_length (int): Seed length to be used in the seed and extend algorithm.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
        exact (bool): If True, will use brute force alignment to align baits instead of seed-and-extend.
        num_cores (int): If greater than 1, will parallelize naive alignment.
    '''
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Calculating coverarage redundancy with provided arguments:\n')
        f.write('d (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(rc))
        f.write('Exact matching = {}\n'.format(exact))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = utils.calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    alignment_function = utils.naive_alignment # for using exact matching
    if not exact:
        alignment_function = utils.seed_and_extend # otherwise
    
    length = len(s)
    cov = [0] * length
    l = len(baits[0])
    ignore = utils.initialize_ignore_vector(seqlens, l)
    if verbose:
        f.write('Initialized integer array and ignore vector with length {}\n'.format(length))

    for id, bait in zip(ids, baits):
        if verbose:
            f.write('Aligning bait {}.\n'.format(id))
        matches = alignment_function(s, bait, d, ignore, sa = sa, num_cores = num_cores)
        coverages = utils.calculate_coverage(matches, l)
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

def coverage_to_workload(id, expected_cov):
    length = len(expected_cov)
    expected_poly = Polynomial([1 - expected_cov[0], expected_cov[0]]) # expected coverage
    for i in range(1, length):
        if expected_cov[i] == 0:
            continue
        curr_poly = Polynomial([1 - expected_cov[i], expected_cov[i]])
        expected_poly = expected_poly * curr_poly
    expected_coeffs = expected_poly.coefficients
    workload = 0
    for i in range(len(expected_coeffs)):
        workload += i * expected_coeffs[i]
    return id, workload

def calculate_workloads(
    ids: List[str],
    baits: List[str],
    s: Union[str, List[str]],
    d: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
    exact: bool = False,
    num_cores = 1
):
    '''
    Given a list of sequences and a solution bait set, returns a dictionary that indicates
    the expected workload of each bait.
    Accuracy of results repend on seed length. For consistency, use the same seed length
    provided to the picking algorithm.

    Arguments:
        ids (list[str]): IDs of solution baits.
        baits (list[str]): Solution baits.
        s (list[str]): Set of sequences that the solution baits cover.
        d (int): Allowed mismatches for bait binding.
        seed_length (int): Seed length to be used in the seed and extend algorithm.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
        exact (bool): If True, will use brute force alignment to align baits instead of seed-and-extend.
        num_cores (int): If greater than 1, will parallelize naive alignment.
    '''
    
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Calculating workloads with provided arguments:\n')
        f.write('d (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(seed_length))
        f.write('Exact matching = {}\n'.format(exact))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = utils.calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    length = len(s)
    l = len(baits[0])
    ignore = utils.initialize_ignore_vector(seqlens, l)

    sa = divsufsort(s)
    
    alignment_function = utils.naive_alignment # for using exact matching
    if not exact:
        alignment_function = utils.seed_and_extend # otherwise

    coverage_vectors = {}
    workloads = {}
    for id, bait in zip(ids, baits):
        coverage_vectors[id] = np.zeros(length)
        if verbose:
            f.write('Aligning bait {}.\n'.format(id))
        matches = alignment_function(s, bait, d, ignore, sa = sa, num_cores = num_cores)
        coverages = utils.calculate_coverage(matches, l)
        if verbose:
            f.write('Bait covers between:\n {}\n'.format(coverages))
        for c in coverages:
            for i in range(c[0], c[1]):
                coverage_vectors[id][i] = 1
    print('Calculated coverages!')
    coverage_sum = np.sum(list(coverage_vectors.values()), axis = 0)
    where = np.where(coverage_sum == 0)[0]
    print('Uncovered indices: ' + str(where))
    for i in where:
        if ignore[i] == True:
            print(str(i) + ' was meant to be ignored.')
    
    coverage_sum[coverage_sum == 0] = 1
    # for id, bait in tqdm(zip(ids, baits), desc = 'Computing workloads', unit = 'baits'):
    args = [(id, coverage_vectors[id] / coverage_sum) for id in ids]
    with Pool(num_cores) as pool:
        results = pool.starmap(coverage_to_workload, args)
    for result in results:
        workloads[result[0]] = result[1]
    return workloads
    
def calculate_coverages(
    ids: List[str],
    baits: List[str],
    s: Union[str, List[str]],
    d: int,
    seed_length: int = 10,
    rc: bool = False,
    log: str = None,
    exact: bool = False,
    num_cores = 1
):
    '''
    Given a list of sequences and a solution bait set, returns a dictionary that indicates
    the coverage of each provided bait.
    Accuracy of results repend on seed length. For consistency, use the same seed length
    provided to the picking algorithm.

    Arguments:
        ids (list[str]): IDs of solution baits.
        baits (list[str]): Solution baits.
        s (list[str]): Set of sequences that the solution baits cover.
        m (int): Allowed mismatches for bait binding.
        seed_length (int): Seed length to be used in the seed and extend algorithm.
        rc (bool): If True, reverse complements will be used for calculating coverage.
        log (str): If not None, the process will log its progress to the specified path.
        exact (bool): If True, will use brute force alignment to align baits instead of seed-and-extend.
        num_cores (int): If greater than 1, will parallelize naive alignment.
    '''
    if log is not None:
        verbose = True
        f = open(log, 'w')
        f.write('Calculating fairness with provided arguments:\n')
        f.write('m (mismatch allowance) = {}\n'.format(d))
        f.write('seed_length = {}\n'.format(seed_length))
        f.write('rc (reverse complements) = {}\n'.format(rc))
        f.write('Exact matching = {}\n'.format(exact))
        f.write('--------\n')
    else:
        verbose = False

    if isinstance(s, list):
        seqlens = utils.calculate_seqlens(s)
        s = ''.join(s)
    else:
        seqlens = [len(s)]

    sa = divsufsort(s)
    
    alignment_function = utils.naive_alignment # for using exact matching
    if not exact:
        alignment_function = utils.seed_and_extend # otherwise
        
    l = len(baits[0])
    ignore = utils.initialize_ignore_vector(seqlens, l)

    result = {}
    for id in ids:
        result[id] = 0

    for i, bait in enumerate(baits):
        if verbose:
            f.write('Aligning bait {}.\n'.format(ids[i]))
        l = len(bait)
        matches = alignment_function(s, bait, d, ignore, sa = sa, num_cores = num_cores)
        coverages = utils.calculate_coverage(matches, l)
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
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        exact = int(sys.argv[7]) > 0
        num_cores = int(sys.argv[8])
        ids, baits = read_sequences(bpath)
        _, sequences = read_sequences(spath)
        vec = calculate_redundancy(
            ids,
            baits,
            sequences,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log',
            exact = exact,
            num_cores = num_cores
        )
        np.savetxt(output, np.array(vec))
        # plt.plot(vec, '-')
        # plt.xlabel('Indices')
        # plt.ylabel('Coverage')
        # plt.yticks(range(max(vec) + 1))
        # plt.title('Redundancy analysis of {} with {}'.format(
        #     os.path.basename(spath), 
        #     os.path.basename(bpath)
        # ))
        # plt.savefig(output)
    elif option == 'workloads':
        bpath = sys.argv[2]
        spath = sys.argv[3]
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        exact = int(sys.argv[7]) > 0
        num_cores = int(sys.argv[8])
        ids, baits = read_sequences(bpath)
        _, sequences = read_sequences(spath)
        workloads = calculate_workloads(
            ids,
            baits,
            sequences,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log',
            exact = exact,
            num_cores = num_cores
        )
        
        # Write the dictionary to a CSV file
        with open(output, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Bait ID', 'Expected Workload'])
            for key, value in workloads.items():
                writer.writerow([key, value])
                
    elif option == 'verify':
        bpath = sys.argv[2]
        spath = sys.argv[3]
        d = int(sys.argv[4])
        output = sys.argv[5]
        rc = int(sys.argv[6]) > 0
        exact = int(sys.argv[7]) > 0
        num_cores = int(sys.argv[8])
        ids, baits = read_sequences(bpath)
        _, sequences = read_sequences(spath)
        uncovered = verify_baits(
            ids,
            baits,
            sequences,
            d,
            seed_length = 10,
            rc = rc,
            log = output[:output.rfind('.')] + '.log',
            exact = exact,
            num_cores = num_cores
        )
        
        # Write the dictionary to a CSV file
        if len(uncovered) > 0:
            print(os.path.basename(bpath), uncovered)
    else:
        print('Invalid option:', option)
    
        
    