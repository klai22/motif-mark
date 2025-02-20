#!/usr/bin/env python

#./Lai_deduper.py -u STL96.txt -f "input_sorted.sam" -o "filtered_sorted.sam"

# Importing modules
import os # works with file names 
import re # a module that provides regex 
import cairo # importing pycairo for vis. 

# MANAGING INPUTS 
#setting get_args 
import argparse
#defining input args 
def get_args():
    parser = argparse.ArgumentParser(description="Motif-mark script to visualize motifs on given sequence(s). This code requires following input: (1) .txt file listing motifs of interest & (2) .fasta file of sequences to search through.")
    parser.add_argument("-f", required=True, help="Absolute path to INPUT fasta file") 
    parser.add_argument("-m", required=True, help="Absolute path to INPUT motifs file")
    return parser.parse_args()

#setting input args--> variables 
args=get_args()
input_fasta_file=args.f
input_motifs_file=args.m

#setting output file name 
input_filename=os.path.basename(input_fasta_file) #setting basename to do splitting on 
filename_prefix,fasta_suffix=os.path.splitext(input_filename) # splitting up base name by "Figure_1.fa" -->("Figure_1",".fa")
output_png_file = f"{filename_prefix}.png" # Setting output filename (adding isolated prefix from previous step + .png)

# READ & PARSE FASTA FILE 
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

# IDENTIFYING MOTIFS in SEQ.s
class Motif: 
    def _init(self, motif_seq):
        "stores a motif seq. + searches for motif in a given seq." 
    
    #Data
    self.motif_seq = motif_seq
    self.positions = [] # empty list to hold positions/coordinates of validated motifs 

    #Methods 
    def find_motif_in_seq(self, seq): 
        "Finds motif in a given sequence, accounts for ambigious nucleotides"
        
        #creating a dictionary to 'translate' ambiguous nucleotides (https://en.wikipedia.org/wiki/Nucleic_acid_notation)
        ambiguous_nucleotides = {
            #purines & pyramadines 
            'R':['A','G'],
            'Y':['C','T'],
            #strong & weak 
            'S':['G','C'],
            'W':['A','T'],
            #amino & ketone 
            'K':['G','T'],
            'M':['A','C'],
            #excluding 1 base 
            'B':['C','G','T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            #any of the 4 bases 
            'N': ['A', 'C', 'G', 'T']
        } 

        # generate all possible motifs (considering ambigious bases)
        def generate_possible_motifs(motif):
            motifs = [motif] #create a permanent list holding each of the possible motif combinations 
            for i,base in enumerate(motif): # loop through each letter in a motif seq. 
                if base in ambiguous_nucleotides: # if there's a base that matches to ambiguous_nucleotides dict....
                    new_motifs = [] #create a temporary list to hold possible motif seq.s 
                    for m in motifs: # for every motif variant that has been generated so far... 
                        for substitute_base in ambiguous_nucleotides[base]: # for every possible substitute base for a given ambiguous base...
                            #adds newly generated motif variant (done below) to new_motifs temporary list. 
                                # write a new string that replaces base @ position i w/ each possible substitute base 
                            new_motifs.append(m[:i] + substitute_base + m[i+1:])
                    motifs = new_motifs # update motifs (permanent list) w/ newly generated motifs 
            return motifs 

        possible_motifs = generate_possible_motifs(self.motif_seq)
    
        #Search for Motifs in Seq. 
        positions = [] # list holding positions/coordinates where matches are detected 
        for motif in possible_motifs:
            for m in re.finditer(motif, seq): # for every combination of motif-var. & seq. 
            # re.finditer returns the position-coordinates where a match was found & stores it in 'm'
                positions.append((m.start(), m.end())) # append coordinates to positions list 
        self.positions = positions # update positions list 
        return positions 

    def find_exons_introns_in_seq(self,seq): 
        "Identifies which regions are exons (uppercase) vs. introns (lowercase) in given seq. Returns a list of tuples w/ (exon/intron, start pos., end pos.) "
        # init variables 
        regions = [] # list to hold output (described in doctstring)
        start = 0 #holds start pos. 
        exon_intron = None # holds strings to differentiate exon vs. intron label assignment 

        for i,base in enumerate(seq): # for every base in sequence 
            #EXON 
            if base.isupper(): #if base is uppercase.....
                # if current region isn't already identified as an exon (exon_intron variable doesn't already have the label 'exon')....
                if exon_intron != 'exon': 
                    # 'Close' previous region if any (add the previously stored (intron) information to regions [] now that we ecountered the start of an exon)
                    if exon_intron: 
                        regions.append((exon_intron,start,i)) 
                    #set exon_intron = 'exon' string now 
                    exon_intron = 'exon' 
                    #reset start coordinate
                    start = i 
            #INTRON 
            elif base.islower(): # if base is lowercase.....
                ## if current region isn't already identified as intron.... 
                if exon_intron != 'intron':
                # 'Close' previous region (exon) if any --> regions [] 
                    if exon_intron: 
                        regions.append((exon_intron,start,i)) 
                    #set exon_intron = 'intron' string now 
                    exon_intron = 'intron'
                    #reset start coordinate 
                    start = i 
        # Encounter / Add final region to regions [] 
        if exon_intron: 
            regions.append((exon_intron,start,len(seq))) 
        
        return regions 


# VISUALIZE MOTIFS IN SEQUENCE 
class Visualizer: 
    def __init__(self,width = 1000, height = 200, scale_factor = 10 ): 
        "Visualizing where the motifs fall amongst the sequence(s) via pycairo"
        #setting parameters set above 
        self.width = width 
        self.height = height 
        self.scale_factor = scale_factor 






