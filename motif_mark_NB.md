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

TO DO: 

*Table of Contents: key tasks* 
* Parse the FASTA file. (extract each seq. --> obj. )
* Parse the motifs file.( search for motif (all its possible versions) in seq. obj.s --> position )
* Handle ambiguous nucleotides.
* Handle overlapping motifs.
* Denote introns and exons.
* Scale features appropriately.
* Generate the final PNG image.

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
2/24/25
____________________________________
# Writing Class to 3. Visualize Motifs on Sequence via Pycairo 

# Write the main() fxn to actually make objects from the classes, and then have them interact to create the visualization 

# Attempting to Execute motif-mark_oop.py 
Command to execute script: 
```
./motif-mark-oop.py -f "/Users/kennethlai/Desktop/BI625/Assignments/OOP/Figure_1.fasta" -m "/Users/kennethlai/Desktop/BI625/Assignments/OOP/Fig_1_motifs.txt"
```
Changing permissions to .py file 
```
chmod 755 motif-mark-oop.py
```
# Debugging Code 
##### Problem 1: 
```
Traceback (most recent call last):
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 31, in <module>
    class Sequence: 
    ^^^^^^^^^^^^^^^
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 35, in Sequence
    self.header = header #FASTA seq. header ">INSR......." 
                  ^^^^^^
NameError: name 'header' is not defined
```
Solution: 
* I forgot to indent the portion below UNDER def __init__(self...) for class Sequence: 
```
self.header = header  # FASTA seq. header ">INSR......."
self.sequence = sequence  # nucleotide seq.
```
##### Problem 2: 
```
Traceback (most recent call last):
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 70, in <module>
    class Motif: 
    ^^^^^^^^^^^^
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 75, in Motif
    self.motif_seq = motif_seq
                     ^^^^^^^^^
NameError: name 'motif_seq' is not defined
```
Solution: 
* same issue as Problem 1 for class Motif: 
```
        #Data
        self.motif_seq = motif_seq
        self.positions = [] # empty list to hold positions/coordinates of validated motifs 
```
##### Problem 3: 
```
Traceback (most recent call last):
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 271, in <module>
    main()
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 239, in main
    sequences = sequence_obj.parse_fasta(input_fasta_file) 
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: Sequence.parse_fasta() takes 1 positional argument but 2 were given
```
Solution: 
* Currently the way .parse_fasta() from class Sequence() is written (below), it only accepts one argument which is the input fasta file. 
```
class Sequence: 
    def __init__(self,header,sequence):
        "stores a fasta sequence"
    #Data
        self.header = header #FASTA seq. header ">INSR......." 
        self.sequence = sequence #nucleotide seq.  

    #Methods 
    def parse_fasta(fasta_file):
        "FASTA file --> list of seq. obj.s"
        sequences = []  #init empty list for each seq. obj. 
        with open(fasta_file,'r') as file:  
            #init. variables to hold headers & seq.s 
            header = None 
            sequence = "" 
            # loop through each line in file 
            for line in file:
                line = line.strip() #strip new line 
                # if encounter a new header line...
                if line.startswith(">"):  
                    # if header is NOT None (collected a seq. from previous loop/line)....
                    if header: 
                        # add (obj. created below) to sequences list (of obj.s) 
                            # create a obj. using Sequence class template to store previous header & sequence 
                        sequences.append(Sequence(header,sequence)) 
                    # set header to current header & reset seq. for next seq. 
                    header = line 
                    sequence = "" 
                # if encounter the sequence line (after header loop)....
                else: 
                    sequence += line # keep adding seq. to current sequence variable/string until hit next header...
                # encountering the last sequence line....
                if header: #if you have a header that isn't = None .... (which will be true for the last line)
                    # add (obj. created below) to sequences list (of obj.s) 
                            # create a obj. using Sequence class template to store previous header & sequence 
                        sequences.append(Sequence(header,sequence)) 
            return sequences # returns list of seq. objects created/parsed from fasta file 
```
* The reason for "TypeError: Sequence.parse_fasta() takes 1 positional argument but 2 were given", is because I forgot to put "self" as an argument in the parse_fasta() method. The udpated version is below 
```
#Methods 
    def parse_fasta(self,fasta_file):
        "FASTA file --> list of seq. obj.s"
        sequences = []  #init empty list for each seq. obj. 
        with open(fasta_file,'r') as file:  
```
##### Problem 4: 
```
Traceback (most recent call last):
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 271, in <module>
    main()
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 265, in main
    visualizer_obj.visualize_seq(scaled_regions,scaled_motifs) # output filename already specified in the class to match the input file's prefix 
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/kennethlai/Desktop/BI625/Assignments/OOP/motif-mark/./motif-mark-oop.py", line 188, in visualize_seq
    surface = cairo.ImageSurface(cairo.FORMAT_RBG24, self.width, self.height)
                                 ^^^^^^^^^^^^^^^^^^
AttributeError: module 'cairo' has no attribute 'FORMAT_RBG24'. Did you mean: 'FORMAT_RGB24'?
```
Solution: 
```
surface = cairo.ImageSurface(cairo.FORMAT_RGB24, self.width, self.height)
```
____________________________________
2/26/25
____________________________________

* Re-wrote code so that Visualizer() outputs are instead saved into a list that holds multiple graphics (PyCairo objects = surfaces) 
* Writing a Combined_Visualization() - type class that can take the graphics list from Visualization() --> paste all into one master .png file. 
    * rewrote Visualizer() & Main() classes to account fo this 

* pushing to GitHub