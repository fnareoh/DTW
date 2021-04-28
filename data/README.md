# Simulation de reads longs: 

**Référence**

Creation d'une séquence de coli de taille 10kb qui servira de référence : 

```bash
fastahack -i /Users/pierre/data/references/ecoli_one_line.fasta
fastahack -r ecoli:100000-110000 /Users/pierre/data/references/ecoli_one_line.fasta > ecoli_10kb.fa 
```

(avec fastahack: `https://github.com/ekg/fastahack`) et ecoli: genome public. 



Puis ajout d'un ligne de header dans le .fa (sinon le format n'est pas respecté)



**Simulation**

Utilisation de l'outil badread: 

```bash
git clone --recursive https://github.com/rrwick/Badread
```



```bash
python Badread/badread-runner.py simulate --length 500,500 --error_model nanopore2020 --start_adapter 0,0  --end_adapter 0,0 --junk_reads 0 --random_reads 0  --chimeras 0 --glitches 0,0,0  --quantity 10x --reference ecoli_10kb.fa | gzip > reads_coli.fastq.gz
```



Le header de chaque read indique entre autres la position de simulation.



J'ai simplifié la simulation ici, mais il est possible d'ajouter des problèmes de séquençage (reads random, portions de séquence random, ...)



**Vérification par remapping sur la référence**

```bash
minimap2 -a ecoli_10kb.fa reads_coli.fastq.gz |gzip > align_reads_coli.sam.gz
```

https://github.com/lh3/minimap2



Le .sam indique pour chaque read comment il a été mappé. 



