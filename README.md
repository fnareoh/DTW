# Frame DTW package

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

Experiments for the default parameters should be reproducible with the following (simplistic) bash script.

```
./experiments.sh
```

# Module description

The source code is available in the `src`folder and organized in sub-modules:
* `dynamic_prog` contains dynamic programming solutions for computing global and local alignments between strings under Edit distance, Homoedit distace, DTW.
* `frameDTW` is a work in progress where the goal is to have a O(knm) time computation of the dynamic programming matrix for DTW distances up to k on a text with n runs and a pattern with m runs.
* `experiments` contains all scripts for comparing how edit distance and DTW distance are affected by homopolymer errors. `read_generator.py` generates reads and computes the ED and DTW distances from them to a reference genome, and `plot.py` compares the resulting distances graphically using matplotlib.

The folder test contain a few testing scripts such as `validation.py` which compares the dynamic programming and frame implementation of local alignment of DTW.
