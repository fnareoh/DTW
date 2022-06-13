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
* `dynamic_prog` contains dynamic programming solution to computing global and local alignment between string.
* `frameDTW` is a work in progress where the goal is to have a O(knm) time computation of the dynamic programming matrix for DTW distances up to k on a text with n runs and a pattern with m runs.
* `experiments` contains every scripts for the actual experiments on the compared efficiency of DTW local alignment vs ED local alignment. In particular, `read_generator.py` generate the reads and compute their local alignment against the genome, where as `plot.py` uses matplotlib to create graphs visualizing those alignment scores.

The folder test contain a few testing scripts such as `validation.py` which compares the dynamic and frame implementation of local alignment of DTW.