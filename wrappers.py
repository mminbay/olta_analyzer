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
    if option == 'vsearch':
        subs_path = sys.argv[2]
        vsearch_id = sys.argv[3]
        basename = os.path.basename(subs_path)[:os.path.basename(subs_path).rfind('.')]
        output_dir = os.path.join(CLUSTER_OUTPUT_DIR, basename)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        command = [
            VSEARCH_PATH,  
            '--cluster_fast',
            subs_path, 
            '--clusters',
            os.path.join(output_dir, 'cluster.fasta'),
            '--id',
            vsearch_id,
            '--gapopen',
            '99999',
            '--gapext',
            '99999',
            '--log',
            os.path.join(CLUSTER_OUTPUT_DIR, basename + '.log')
        ]
        print('Starting to cluster substrings...')
        start = time.perf_counter()

        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)

        process.wait()
        end = time.perf_counter()
        print('Entire process took {} seconds'.format(str(end - start)))
    elif option == 'catch':
        seq_path = sys.argv[2]
        bait_length = sys.argv[3]
        mismatch = sys.argv[4]
        basename = os.path.basename(seq_path)[:os.path.basename(seq_path).rfind('.')]
        output_dir = os.path.join(SOFTWARE_OUTPUT_DIR, basename)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        command = [
            'python',
            'C:\\Users\\mminbay\\AppData\\Local\\anaconda3\\envs\\bait_research\\Scripts\\design.py',  
            seq_path,
            '-o', 
            os.path.join(output_dir, basename + '_catch.fasta'),
            '-pl',
            bait_length,
            '-m',
            mismatch,
        ]
        print('Starting to pick baits with CATCH...')
        start = time.perf_counter()

        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)

        process.wait()
        end = time.perf_counter()
        print('Entire process took {} seconds'.format(str(end - start)))
    elif option == 'synthesize':
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
        