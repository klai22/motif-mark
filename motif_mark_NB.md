Kenny's Notebook for OOP Motif Mark Assignment (BI625) 

Assignment: https://github.com/klai22/motif-mark.git

Data / Notebook/ WD Location: /Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark

Creating interactive sessions: srun --account=bgmp -p bgmp -N 1 -c 4 --mem=100G -t 6:00:00 --pty bash

SLURM TRACKER Template: 
| Sample | JobID | Run Time(mm:ss) | CPU Usage (%)| Exit Status |
|---|---|---|---|---|
|6_2D|15903418|1:13.24|671|0|
____________________________________
2/19/25
____________________________________
# Setting up Project. 
### Create a github repo. for assignment: 

```
wd: /Users/kennethlai/Desktop/BI625/Assignments/OOP/
git init 
git clone https://github.com/klai22/motif-mark.git
```
### Set up Conda Env. for pycairo 
```
conda create -n pycairo_env pycairo 
conda activate pycairo_env
conda install pycairo 
```
### Create motif-mark-oop.py
```
nano motif-mark-oop.py
```

# Writing Classes to 1. Parse Fasta file & 2. Parse motifs file. 
Wrote classes to hold Fasta info. & Motifs info. Left off on p.3 - a visualizer class to actually visualize where motifs fall on sequence(s). 

____________________________________
2/____/25
____________________________________
# Writing Class to 3. Visualize Motifs on Sequence via Pycairo 

# --> Integrate exon & intro info. into main code...



TO DO: 

*Table of Contents: key tasks* 
* Parse the FASTA file. (extract each seq. --> obj. )
* Parse the motifs file.( search for motif (all its possible versions) in seq. obj.s --> position )
* Handle ambiguous nucleotides.
* Handle overlapping motifs.
* Denote introns and exons.
* Scale features appropriately.
* Generate the final PNG image.
