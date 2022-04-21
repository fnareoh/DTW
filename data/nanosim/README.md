https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11

Shorten Genome:
```
# To select the first line for Fasta format
head -n 1 NC_000001_11.fasta > short_NC_000001_11.fasta
# To extract the lines 500 to 1500
head -1500 NC_000001_11.fasta | tail -300 >> short_NC_000001_11.fasta
```

Nanosim must be installed with conda before
```
conda activate nanosim
#Fasta output
./src/simulator.py genome -rg short_NC_000001_11.fasta -c pre-trained_models/human_NA12878_DNA_FAB49712_albacore/training -n 200
#Fastq output
./src/simulator.py genome -rg short_NC_000001_11.fasta -c pre-trained_models/human_NA12878_DNA_FAB49712_albacore/training -n 200 --fastq --basecaller albacore --max_len 500
```

Minimap2 alignement
```
minimap2 -a short_NC_000001_11.fasta short_simulated_aligned_reads.fastq  >short_simulated_aligned_reads.sam
```