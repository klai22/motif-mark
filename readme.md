Kenneth Lai 
BI625: OOP Motif Mark 
3/6/25 

# Overview 
```motif-mark-oop.py``` is a Python-based algorithim that visualizes motifs on sequences. This tool requires a FASTA file input (holding multiple sequences) and a .txt file listing out the motifs of interest. ```Motif-mark-oop.py``` outputs a master figure (.png format) displaying where the various motifs for this assignment fall on the sequence. The visualization also shows where exons (lower case letters = blue boxes) and introns (upper case letters = black lines) lie on the provided sequence. Furthermore, ambigious nucleotide encoding is supported by this algorithm. 

### Installation: 
* Ensure your have a Python 3 -compatible envrionemnt. Conda package management is recommended. 
1. Create a conda env. with neccesary tools. 
```
conda create -n motif-mark
conda activate motif-mark 
```
2. Clone this repository 
```
git clone https://github.com/klai22/motif-mark/tree/main
cd motif-mark
```

### Executing Algorithm at Command-Line Level 
3. Run the ```motif-mark-oop.py``` script, using ```-f``` & ```-m``` to specify file inputs. 
```
./motif-mark-oop.py -f "/path/to/sequences_file.fasta" -m "/path/to/motifs_file.txt"
```
##### DATA INPUT REQUIREMENTS 
* ```-f```: Input Fasta File: ```sequences_file.fasta`` - holds sequences, each denoted w/ ```>``` as the sequence header. 
* ```-m```: Input Motifs File: ```motifs_file.txt``` - holds motif sequences, separated by ```/n``` s. 

##### DATA OUTPUT 
* Whatever prefix your ```.fasta``` file is labled with will carry over into the output ```.png``` holding the motif visualization. 
* Ex: If your .fasta file is labeled ```Figure_1.fasta```, the output will be labeled ```Figure_1.png```. 


### Goals
* Write a python script utilizing OOP (object-oriented programming) to visualize motifs on multiple sequences. 
* Algorithm can handle FASTA input with sequences that are less than < 1000 bases long & motifs that are < 10 bases long. 
* Include motif, exon, and intron features in graphic while accounting for accurate scaling. Output visualization as a .png file. 
* Handle overlapping motifs. 


### Developer Notes 


* The following encoded was utilized to decipher ambigious nucelotides. 

| Ambigious Nucleotide | Nucleotide Encodings |
|---|---|
Purines & Pyramadines
|R|A,G|
|Y|C,T|
Strong & Weak 
|S|G,C|
|W|A,T|
Amino & Ketone
|K|G,T|
|M|AC|
Excluding 1 base
|B|G,T,C|
|D|A,G,T|
|H|A,T,C|
|V|A,G,C|
Any of the 4 bases
|N|A,G,T,C|

* The key for the visualizations are hard-coded. Thus, only the motifs below are accounted for in the legend(s). 
```
ygcy
GCAUG
catag
YYYYYYYYYY
```
