# Activate the virtual environment
source .venv/bin/activate

# Creates result folder if it doesn't exist already
mkdir -p results;

# Generate reads and compute distances
N=600
# N reads of length 500 for each homopolymer probability
python src/experiments/read_generator.py data/ecoli_10kb.fa -N $N

# Plot scripts
python src/experiments/plot.py results/ecoli_10kb_N_$N ID 0.05