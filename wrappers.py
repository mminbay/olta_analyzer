from fasta_utils import synthesize_sequence, sequences_to_fasta
import os
import subprocess
import sys
import time

VSEARCH_PATH = 'C:\\Users\\mminbay\\Documents\\bait_research\\vsearch-2.24.0-win-x86_64\\bin\\vsearch'
CLUSTER_OUTPUT_DIR = 'C:\\Users\\mminbay\\Documents\\bait_research\\clusters\\'
SOFTWARE_OUTPUT_DIR = 'C:\\Users\\mminbay\\Documents\\bait_research\\outputs\\'

if __name__ == '__main__':
    option = sys.argv[1]
    if option == 'synthesize':
        s_length = int(sys.argv[2])
        n_repeats = int(sys.argv[3])
        repeat_length = [int(x) for x in sys.argv[4].split('-')]
        repeat_coverage = float(sys.argv[5])
        repeat_noise = int(sys.argv[6])
        output_path = sys.argv[7]
        seq = synthesize_sequence(
            s_length,
            n_repeats,
            repeat_length,
            repeat_coverage,
            repeat_noise,
            log = output_path[:output_path.rfind('.')] + '.log' # might want to change this later
        )
        sequences_to_fasta(seq, output_path)
        