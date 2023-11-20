import os
from Bio import SeqIO
import math
from typing import Dict, List
from tqdm import tqdm

class ClusterHelper():
    def __init__(
        self, 
        l: int, 
        k: int,
        seqlens: List[int],
        cluster_map: Dict[int, List[int]]= None,
        cluster_dir: str = None
    ):
        '''
            Arguments:
                seqlens: a list of integers representing the sequence lengths of the sequences to be covered. i'th entry should be length of i'th sequence
                cluster_map: a map that maps representatives to the list of substrings they cover.
                l: length of substrings
                k: used step size when picking substrings to cluster
        '''
        if cluster_map is None and cluster_dir is None:
            raise Exception('A cluster map or directory must be provided')
        if cluster_map is not None and cluster_dir is not None:
            raise Exception('Both a map and directory were provided: Provide one only.')

        if cluster_map is None:
            cluster_map = map_from_clusters(cluster_dir)

        self.seqlens = seqlens
        self.bait_to_subs = cluster_map
        self.sub_to_baits = {}
        last = 0
        for key in cluster_map.keys():
            items = cluster_map[key]
            for item in items:
                self.sub_to_baits[item] = key
        self.lasts = [seqlen - l for seqlen in seqlens]
        self.l = l
        self.k = k

    def get_subs(self, bait: int):
        '''
        Given a substring index of a bait, return the indices
        of the substrings that are covered by that bait.

        Arguments:
            bait (int): Substring index of a bait.
        
        Returns:
            list(int): Indices of substrings covered by that bait.
        '''
        return self.bait_to_subs[bait]

    def get_repr(self, index):
        '''
        Given an index, return the substring indices of the
        baits that cover that index.

        Arguments:
            index (int): An index.
        
        Returns:
            list(int): Substring indices of baits that cover that index.
        '''
        index_seq = 0
        offset = 0
        curr_total_seqlen = self.seqlens[index_seq]
        while index >= curr_total_seqlen:
            offset = curr_total_seqlen
            index_seq += 1
            curr_total_seqlen += self.seqlens[index_seq]

        # calculate the original index, which is the index of this character in its RESPECTIVE sequence (not the entire, concatenated sequence)
        index = index - offset
        # calculate the original index of the last substring that covered this index
        end = (index // self.k) * self.k 
        # calculate the original index of the last substring that covered this index
        start = math.ceil((index - self.l + 1) / self.k) * self.k
        if start < 0:
            start = 0
        
        subs = list(range(start, end + 1, self.k))
        
        result = []
        for sub in subs:
            # add the offset back so we get the concatenated index (instead of the true index)
            sub = sub + offset 
            if sub not in self.sub_to_baits.keys():
                break
            else:
                res = self.sub_to_baits[sub]
                if res not in result:
                    result.append(res)
        if (index + offset) - self.l + 1 >= self.lasts[index_seq] + offset:
            res = self.sub_to_baits[self.lasts[index_seq] + offset]
            if res not in result:
                result.append(res)
        return result

    def get_exact_repr(self, index):
        return self.sub_to_baits[index]
                
def map_from_clusters(dir: str):
    '''
    Reads the output directory of a VSEARCH clustering run, and returns
    a map that maps every cluster representative to the substring indices
    that are represented by it.

    Arguments:
        dir (str): Path to the directory of VSEARCH output.
    
    Returns:
        center_to_cluster (dict(int: int)): A map that maps representative 
            indices to the list of substring indices they represent.
    '''
    center_to_cluster = {}
    
    data = os.listdir(dir)
    for filename in tqdm(data, desc = 'Mapping clusters', unit = 'files'):
        if 'fasta' in filename:
            file_path = os.path.join(dir, filename)

            current_cluster = []
            
            for record in SeqIO.parse(file_path, 'fasta'):
                sequence_name = int(record.id)
    
                current_cluster.append(sequence_name)
    
                if len(current_cluster) == 1:
                    center_sequence_name = sequence_name

            center_to_cluster[center_sequence_name] = current_cluster

    return center_to_cluster