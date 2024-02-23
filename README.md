# bait_research
Work I am doing for my CS honors thesis at colgate.

# Usage
`pick_baits.py` contains all of the methods you will need for generating baits. The file must be run as follows:

```
$ python pick_baits.py <syotti|syotti_s> <spath> <l> <m> <output> <rc>
```

* `syotti`:  Mock implementation of the Syotti heuristic \[1\]. Note that this is a Python implementation with no parallelization - it runs much slower than the original Syotti heuristic and was only implemented for testing and comparison purposes.
* `syotti_s`:  Our algorithm. This expands on the single-pass Syotti heuristic by computing the coverage for every l-length substring that includes the current uncovered index, and picking the substring with the greatest coverage as a bait. This consistently results in fewer baits with a running time no worse than that of **MOCK SYOTTI** multiplied by the bait length. It is in future work to implement this more efficiently to have a running time comparable to the original Syotti implementation.
* `spath`: Path to a `.fasta` file containing the sequences you wish to work with.
* `l`: Positive integer indicating the length of baits to be generated.
* `m`: Non-negative integer indicating the maximum Hamming distance that is permitted for a bait to be considered as "covering" a region.
* `output`: Output path.
* `rc`: Non-negative integer indicating whether reverse complements should be used for calculating bait coverages. Pass 0 to ignore reverse complements, or any number greater to include them.


# Citations
\[1\] J. N. Alanko, I. B. Slizovskiy, D. Lokshtanov, T. Gagie, N. R. Noyes, and C. Boucher, “Syotti: scalable bait design for DNA enrichment,” Bioinformatics, vol. 38, no. Supplement_1, pp. i177–i184, Jun. 2022, doi: 10.1093/bioinformatics/btac226.

