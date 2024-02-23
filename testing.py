import subprocess
import os
import numpy as np
import re
from multiprocessing import Pool

SEQ_DIR = '/home/mminbay/honors_thesis/sequence_data/'
COLI3_SEQ_DIR = os.path.join(SEQ_DIR, 'coli3')
MEGARES_SEQ_DIR = os.path.join(SEQ_DIR, 'megares')
SYNTH_SEQ_DIR = os.path.join(SEQ_DIR, 'synthesized')

OUT_DIR = '/home/mminbay/honors_thesis/outputs'
COLI3_OUR_DIR = os.path.join(OUT_DIR, 'coli3')
MEGARES_OUT_DIR = os.path.join(OUT_DIR, 'megares')
SYNTH_OUT_DIR = os.path.join(OUT_DIR, 'synthesized')

def get_test_info(str):
    def __lookahead(field):
        return '(?=' + re.escape(field) + ')'
    fields = ['L', 'unique', 'l', 'pccov', 'noise']
    result = {}
    for field in fields:
        regex = r'(^|(?<=_))(\d|-)+' + __lookahead(field)
        result[field] = re.search(regex, str)[0]
    return result

def run_test(
    test_path,
    method,
    output_dir,
    d = 40,
    l = 120,
    rc = 0
):
    if method not in ['syotti', 'syotti_smart', 'syotti_wfc', 'wfc_iter']:
        raise Exception('Invalid method')

    file_name = os.path.basename(test_path)
    output_path = os.path.join(output_dir, file_name.replace('.fasta', '_' + method + '.fasta'))
    command = [
        'python', 
        '/home/mminbay/honors_thesis/bait_research/pick_baits.py',
        method,
        test_path,
        str(l),
        str(d),
        output_path,
        str(rc)
    ]
    subprocess.run(command)

def synth_tests(method):
    PROCESSES = 31
    synth_seqs = [file for file in os.listdir(SYNTH_SEQ_DIR) if '.fasta' in file]
    synth_1m = [file for file in synth_seqs if int(get_test_info(file)['L']) <= 1000000]
    synth_1m = sorted(synth_1m, key = lambda x: int(get_test_info(x)['L']))
    synth_1m_paths = [os.path.join(SYNTH_SEQ_DIR, path) for path in synth_1m]

    args = [(
        path,
        method,
        SYNTH_OUT_DIR,
    ) for path in synth_1m_paths]

    with Pool(PROCESSES) as pool:
        result = pool.starmap(run_test, args)

def megares_tests(method):
    PROCESSES = 10
    megares_seqs = [file for file in os.listdir(MEGARES_SEQ_DIR) if '.fasta' in file]
    megares_1m = [file for file in megares_seqs if int(re.search(r'\d+', file)[0]) <= 1000000]
    megares_1m = [os.path.join(MEGARES_SEQ_DIR, path) for path in megares_1m]

    args = [(
        path,
        method,
        MEGARES_OUT_DIR,
    ) for path in megares_1m]

    with Pool(PROCESSES) as pool:
        result = pool.starmap(run_test, args)

def coli3_tests(method):
    PROCESSES = 10
    coli3_seqs = [file for file in os.listdir(COLI3_SEQ_DIR) if '.fasta' in file]
    coli3_1m = [file for file in coli3_seqs if int(re.search(r'\d+', file)[0]) <= 1000000]
    coli3_1m = [os.path.join(COLI3_SEQ_DIR, path) for path in coli3_1m]

    args = [(
        path,
        method,
        COLI3_OUR_DIR,
    ) for path in coli3_1m]

    with Pool(PROCESSES) as pool:
        result = pool.starmap(run_test, args)
    
def main():
    # run_test(
    #     '/home/mminbay/honors_thesis/sequence_data/synthesized/100000L_100unique_120-240l_80pccov_40noise_no5.fasta',
    #     'wfc_iter',
    #     SYNTH_OUT_DIR
    # )
    coli3_tests('wfc_iter')
    coli3_tests('syotti_wfc')
    
if __name__ == '__main__':
    main()