# BlockDTW package

## Install

We first recommend creating a virtual environment:
```
python3 -m venv .venv
```
Then, to activate it, use the following command:
```
source .venv/bin/activate
```

For ease of use the project has been organized as a pip package, for a local editable install use the following command:

```
pip install -e .
```

# Experiment reproduction

Experiments for the default parameters should be reproducible with the following (simplistic) bash script. The results are outputted progressively to a `results` folder.

```
./experiments.sh
```

# Module description

The source code is available in the `src`folder and organized in sub-modules:
* `dynamic_prog` contains dynamic programming solutions for computing global alignments and pattern matching between strings under the edit and DTW distances.
* `BlockDTW` is an implementation of an algorithm that computes the distances between a pattern P of length M with m runs and a text T of length N with n runs in time O(nM+mN), with some extra optimisations for the small-distance regime.  
* `experiments` is a module for comparing how the edit distance and DTW distance are affected by homopolymer errors. In particular, `read_generator.py` generates reads and computes the edit and DTW distances from them to a reference genome, and `plot.py` compares the resulting distances graphically using matplotlib.

The folder test contain a few testing scripts such as `validation.py` which compares the dynamic programming and block implementation of pattern matching for DTW.
