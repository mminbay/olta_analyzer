import sys
from tqdm import tqdm
from cluster_tools import ClusterHelper
from utils import hamming_distance
from fasta_utils import sequences_from_fasta

def stupid_verify_baits(
    seq_path: str,
    baits_path: str, 
    mismatch: int,
):
    '''
    '''
    _, seqs = sequences_from_fasta(seq_path)
    _, baits = sequences_from_fasta(baits_path)

    sequence = ''.join(seqs)
    length = len(sequence)
    l = len(baits[0])
    cov = [False] * length

    data = range(length - l + 1)
    for i in tqdm(data, desc = 'Covering indices', unit = 'item'):
        sub = sequence[i:i + l]
        for bait in baits:
            if hamming_distance(sub, bait) <= mismatch:
                for j in range(l):
                    cov[i + j] = True
                break
    return cov.count(False)

def compare_baits(
    seq_path: str,
    other_baits_path: str,
    cluster_baits_path: str,
    cluster_dir: str,
    l: int,
    k: int,
    mismatch: int,
    output: str,
):
    _, sequences = sequences_from_fasta(seq_path)
    seqlens = [len(seq) for seq in sequences]
    s = ''.join(sequences)
    s = s.upper()
    o_ids, o_baits = sequences_from_fasta(other_baits_path)
    o_baits = [x.upper() for x in o_baits]

    c_ids, c_baits = sequences_from_fasta(cluster_baits_path)
    c_baits = [x.upper() for x in c_baits]
    c_ids = [int(x) for x in c_ids]

    o_baits = list(zip(o_ids, o_baits))
    c_bait_to_index = dict(list(zip(c_baits, c_ids)))

    cluster = ClusterHelper(
        l,
        k,
        seqlens,
        cluster_dir = cluster_dir
    )

    with open(output, 'a') as file:
        file.write(
            'Comparing the baits at {} with the representatives at {}.\n'
            .format(other_baits_path, cluster_baits_path)
        )
        file.write(
            'Using VSEARCH clustering output at {}.\n\n'.format(cluster_dir)
        )
        file.write('Number of baits on comparison file: {}\n'.format(len(o_baits)))
        file.write('Number of representatives on cluster results file: {}\n'.format(len(c_baits)))
        file.write('Difference: {}\n\n'.format(len(c_baits) - len(o_baits)))

        match = 0
        cov_diff = 0
        data = o_baits
        for bait in tqdm(data, desc = 'Comparing baits', unit = 'baits'):
            id = bait[0]
            bait_seq = bait[1]
            file.write('Identifier {}, sequence: {}\n'.format(id, bait_seq))
            if bait_seq in c_bait_to_index.keys():
                match += 1
                file.write('Bait is present in our set of representatives.\n')
                i = c_bait_to_index[bait_seq]
                repr = cluster.get_exact_repr(i)
                file.write('The index of the representative is {}.\n'.format(i))
            else:
                file.write('Bait is NOT present in our set of baits.\n')
                i = s.find(bait_seq)
                repr = cluster.get_exact_repr(i)
                file.write('An exact occurence of this substring is at {}. '.format(i))
                file.write('The index of this occurence\'s representative is {}.\n'.format(repr))
            coverage = cluster.get_subs(repr)
            file.write('Clustering-computed coverage of this representative is {} substrings.\n'.format(len(coverage)))
            file.write('These substrings\'s indices are: \n{}\n'.format(coverage))
            stupid_coverage = []
            for k in range(len(s) - l + 1):
                sub = s[k:k + l]
                if hamming_distance(bait_seq, sub) <= mismatch:
                    stupid_coverage.append(k)
            file.write('Manual alignment coverage of this representative is {} substrings.\n'.format(len(stupid_coverage)))
            file.write('These substrings\'s indices are: \n{}\n\n'.format(stupid_coverage))

            cov_diff += len(stupid_coverage) - len(coverage)

        file.write('--- SUMMARY --- \n')
        file.write('{} baits were present in both lists.\n'.format(match))
        file.write('Cluster-computed coverage missed {} alignments in total'.format(cov_diff))
        

if __name__ == '__main__':
    method = sys.argv[1]
    if method == 'stupid_verify':
        seq_path = sys.argv[2]
        baits_path = sys.argv[3]
        mismatch = int(sys.argv[4])
        print(stupid_verify_baits(seq_path, baits_path, mismatch))
    if method == 'compare':
        seq_path = sys.argv[2]
        o_baits_path = sys.argv[3]
        c_baits_path = sys.argv[4]
        cluster_dir = sys.argv[5]
        l = int(sys.argv[6])
        k = int(sys.argv[7])
        m = int(sys.argv[8])
        output = sys.argv[9]
        compare_baits(
            seq_path,
            o_baits_path,
            c_baits_path,
            cluster_dir,
            l,
            k,
            m,
            output
        )