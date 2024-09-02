# bait_research
Analyzer used for OLTA (https://github.com/FuelTheBurn/generative-bait-clustering)

# Usage
## Getting started
We recommend using a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) to manage the dependencies of the analyzer. Once conda is set up, you can create an environment using the command `conda env create --name NAME --file dependencies.yaml`, where `NAME` is replaced by your preferred name for the environment.

## Postprocessing analyses
`postprocessing.py` has the functions for redundancy and workload analyses. You can run the analyses as follows:

```
$ python postprocessing.py [redundancy, workloads] <bait_path> <seq_path> <mismatch> <output> <rcomp> <exact> <n_cores>
```
Arguments:
* `bait_path`: The path to the baits that are to be analyzed. This should be a `.fasta` or `.fna` file where each entry is a bait.
* `seq_path`: The path to the sequence(s) that the baits will be analyzed on. This should be a `.fasta` or `.fna` file where each entry is an input sequence.
*  `mismatch`: Number of mismatches that are tolerated for bait hybridizations (the $\theta$ parameter).
*  `output`: The path where the analysis result will be outputted. Redundancy analyses require a `.txt` extension file, as the output will be a vector where each entry indicates how many baits in the solution set cover that respective position in the input sequence(s). Workload analyses require a `.csv` extension, as the output will be a table indicating the expected workload of each bait. See the original paper for the definitions of these terms.
*  `rcomp`: Currently non-functional. An integer greater than 0 will consider reverse complement alignments of the baits.
*  `exact`: Using an integer greater than 0 will use the brute-force algorithm for calculating bait alignments. Using 0 will use the seed-and-extend alignment. We recommend using the brute-force algorithm when analyzing baits produced by OLTA.
*  `n_cores`: Number of processes to be used if using the naive alignment.

### Example

`example_bait` contains a bait set that was produced by OLTA for the synthetic sequence `synth_data/L500000_RL(120,240)_RC50_RE40_no0.fasta`. To do a redundancy analysis with exact alignment using 8 processors, we can run the following command:
```
$ python postprocessing.py redundancy example_bait/L500000_RL\(120\,240\)_RC50_RE40_no0_baits.txt synth_data/L500000_RL\(120\,240\)_RC50_RE40_no0.fasta 40 example_output.txt 0 1 8
```
The resulting numpy array will be outputted to `example_output.txt`.

## Sequence synthesizing
`wrappers.py` wraps the function for synthesizing sequences with controlled repetitions. You can use this functionality as follows:
```
$ python wrappers.py synthesize <sl> <rn> <rl> <rc> <modif> <output>
```
Arguments:
* `sl`: Length of the sequence to be synthesized in nucleotides.
* `rn`: Number of unique seed repeats to be planted in the sequence.
* `rl`: Minimum and maximum possible lengths (in nucleotides) for the unique seed repeats, separated by a hyphen (-).
* `rc`: Fraction of the sequence to be covered with imperfect copies of the unique seed repeats.
* `modif`: Number of modifications to be made on a seed repeat before planting it into the sequence.
* `output`: The path where the resulting sequence will be outputted.

### Example
```
python wrappers.py synthesize 50000 50 120-240 0.5 40 test.fasta
```
Running this command will generate a sequence with 50000 nucleotides. Half of this sequence will be covered by imperfect copies of 50 seed repeats, each of which are 120-240 nucleotides in length. The imperfect copies will have at most 40 Hamming distance to their respective seeds.

### Notes
As the repeat coverage grows beyond 0.5, it becomes increasingly likely that the algorithm will result in a dead end, i.e. there will be no more contiguous empty regions to place repetitions despite not having reached the desired fraction. Multiple attempts might be necessary to synthesize such sequences.

Using a repeat coverage of 1.0 will just concatenate the imperfect copies instead of planting them into an empty sequence.