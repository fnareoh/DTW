# Creation of ecoli_10kb.fa

Creation of an E.coli of size 10kb which will be used as a reference Genome from which we generate a "mutated" genome and some reads : 

```bash
fastahack -i /Users/pierre/data/references/ecoli_one_line.fasta
fastahack -r ecoli:100000-110000 /Users/pierre/data/references/ecoli_one_line.fasta > ecoli_10kb.fa
```



