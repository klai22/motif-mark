#!/usr/bin/env python

#./motif-mark-oop.py -f "/Users/kennethlai/Desktop/BI625/Assignments/OOP/Figure_1.fasta" -m "/Users/kennethlai/Desktop/BI625/Assignments/OOP/Fig_1_motifs.txt"

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

# CLASS: READ & PARSE FASTA FILE 
class Sequence: 
    def __init__(self,header,sequence):
        "stores a fasta sequence"
    #Data
        self.header = header #FASTA seq. header ">INSR......." 
        self.sequence = sequence #nucleotide seq.  

    #Methods 
    def parse_fasta(self,fasta_file):
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
    # Defining what shows up when you want to print the object (for DEBUGGING)
    def __str__(self):
        return f"Header: {self.header}\nSequence: {self.sequence[:50]}..."  # Print first 50 char.

# CLASS: IDENTIFYING MOTIFS in SEQ.s
class Motif: 
    def __init__(self, motif_seq):
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
                positions.append((m.start(), m.end(),self.motif_seq)) # append coordinates to positions list 
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


# CLASS: VISUALIZE MOTIFS IN SEQUENCE 
class Visualizer: 
    def __init__(self,width = 2900, height = 1100, scale_factor = 3.5, left_margin_offset=20,bottom_margin_offset=90): 
        "Visualizing where the motifs fall amongst the sequence(s) via pycairo"
        #setting parameters set above 
        self.width = width 
        self.height = height 
        self.scale_factor = scale_factor 
        self.left_margin_offset = left_margin_offset
        self.bottom_margin_offset = bottom_margin_offset

    def scaling_feature_coordinates(self, regions, motif_positions):
        "Multiplies & maps start & end feature-coordinates (regions [exon & introns] + motifs) to scale_factor # so that the image will fit within self's image's specified width & height. Returns same updated list-objects, but with re-scaled position numbers."
        #scaling region coordinates
            # applying the following scaling-conversions (multiplications) (exon_intron, start * self.scale_factor, end * self.scale_factor) for every tuple in the regions [] 
            ## adding a left margin so that image isn't touching left side of image
        scaled_regions = [(exon_intron, (start + self.left_margin_offset) * self.scale_factor, (end + self.left_margin_offset) * self.scale_factor)
                          for exon_intron, start, end in regions]
        #scaled_regions = [(exon_intron, start * self.scale_factor, end * self.scale_factor)for exon_intron, start, end in regions]
        #scaling motif coordinates
            # applying the following scaling-conversions (multiplications) (start * self.scale_factor, end * self.scale_factor) for every tuple in the motif_positions [] 
            # adding a left margin so that image isn't touching left side of image
        scaled_motif_positions = [((start + self.left_margin_offset) * self.scale_factor, (end + self.left_margin_offset) * self.scale_factor,motif_seq)
                          for start, end, motif_seq in motif_positions]
        #scaled_motif_positions = [(start * self.scale_factor, end * self.scale_factor)for start, end in motif_positions]
        return scaled_regions, scaled_motif_positions

    def visualize_seq(self, scaled_regions, scaled_motif_positions, seq_title):
        "Visualizes sequence (exons, introns, motifs) via PyCairo. Returns output surface (pyCairo object)"
        #create the coordinates/dimensions to display your graphic, designate output
            #RGB24 specifies RGB format --> 24 bits per pixel (8 bits for each the 3 color channels: Red, Green, Blue).
        surface = cairo.ImageSurface(cairo.FORMAT_RGB24, self.width, self.height)
        #create the coordinates you will be drawing on (like a transparency) - you can create a transformation matrix --> think of this as the "pen" used to add shapes, lines, & text to image
        context = cairo.Context(surface) 
        #set white background 
        context.set_source_rgb(1,1,1) 
        context.paint() #fill entire image surface w/ the color specified above (white)

        # Visualization Parameters 
            # Adjusting y-position of motifs to appear above sequence line & so sequence line appears at very bottom
        sequence_line_y = self.height - self.bottom_margin_offset # Place Sequence line near bottom of image 
        motif_height = 22 # height of motif boxes 
        motif_margin = 10 # distance b/w 1st motif layer & sequence (exons & introns) line 

        # Sorting (ordering) scaled_regions by exon vs. intron (x[0]) so that all exon boxes cover the intron line. 
        #scaled_regions = sorted(scaled_regions,key=lambda x:x[0])

        # Draw Exon & Intron Regions 
        # (A) Draw Introns First (so that they don't cover exon boxes )
        for region in scaled_regions: # loop through exon_intron region coordinates
            exon_intron, start, end = region # separate out the given region tuple 
            # INTRONS = draw a black line (smaller rectangle) 
            if exon_intron == 'intron':
                context.set_source_rgb(0,0,0) # black color 
                context.move_to(start,sequence_line_y+26) # Start point [y = sequence_line_y+26, so midpoint lines up w/ exon's ]
                context.line_to(end,sequence_line_y+26) # End point [y = sequence_line_y+26, so midpoint lines up w/ exon's ]
                context.set_line_width(2) # line width = 2 
                context.stroke() # draws line  
        # (B) Draw Exons Second (so that they don't get covered by intron lines )
        for region in scaled_regions: # loop through exon_intron region coordinates
            exon_intron, start, end = region # separate out the given region tuple 
            # EXONS = draw a blue rectangle 
            if exon_intron == 'exon':
                # set EXON color = (light) BLUE  
                context.set_source_rgb(0.25, 0.60, 0.82)
                # set EXON shape = rectangle & specifying its coordinates 
                context.rectangle(start, sequence_line_y , end - start, 50) 
                    #start & end = scaled position coordinates for exon (x left-side)
                    #sequence_line_y = y top-side 
                    #end - start = width of rectangle 
                    #sequence_line_y= height of rectangle
                # fill rectangle w/ blue color 
                context.fill() 

        # Draw Motifs Positions + slightly offset each one to account to handle overlaps 
        # CREATING COLOR (RGB) MAPPINGS (a dict.)
        motif_color_map = {
            '___________':(1,1,1), # white
            'REGIONS':(1,1,1), # white
            'Exon':(0.25, 0.60, 0.82), # light blue 
            'Intron':(0,0,0), #black
            '__________':(1,1,1), # white
            'MOTIFS':(1,1,1), # white
            'ygcy':(141/255, 48/255, 217/255), # purple
            'GCAUG':(217/255, 144/255, 48/255), # orange
            'catag':(217/255, 48/255, 163/255), # pink
            'YYYYYYYYYY':(53/255, 184/255, 88/255), # green
            '____________':(1,1,1) # white
            }

        # VISUALIZE MOTIFS WHILE ACCOUNTING FOR OVERLAP 
        # 0. Init. list to track ongoing motif layers (visualization of overlapping motifs will be actively tracked)
        motif_layers = [] #stores (end position, layer index) for each motif that is actively overlapping with the current motif being processed. 

        # 1. Sort motifs by Start Position to process them in order 
            #sorted() sorts the list of tuples in INCREASING order 
            #key specifies sorting criteria 
                #lambda x:x[0] = returns first value (start position) of a given tuple (x) 
        sorted_motifs = sorted(scaled_motif_positions, key = lambda x:x[0])
        # 2. Loop through sorted motifs 
        for motif in sorted_motifs: # loop through SORTED motif_position coordinates 
            start,end,motif_seq = motif # separate out the given motif coordinates 
        # 3. Remove motifs that no longer overlap ( if there are any atm) / only keep motifs that overlap (end positions are > current motif's start position )
            motif_layers = [layer for layer in motif_layers if layer[0] > start] # loop through each layer in motif_layers & isolate/keep only motif (tuples) whose end position (layer[0]) is > than the start position of the current motif! 
            # Filtering out motifs that do not overlap frees up layers (space) for other motifs that do not overlap with the current one's end position. 
        # 4. Assign the lowest-possible / next-available layer (y-position) for the given motif. 
            layer = len(motif_layers) # this makes sense because if there are 3 overlapping motifs at the given moment, you would want to assign it a y-coordinate that places it in the 3rd layer so that it sits above the previous 2 already overlapping motifs.
            motif_layers.append((end,layer)) # append motif's end position & assigned layer # to motif_layers [] 
        # 5. Offset motif's y-coordinate to avoid overlap 
            offset = layer * 35 # offset is unique to layer so that each motif coordinate is uniquely offset according to how many overlaps occuring for that given motif atm 
        # 6. MOTIF = draw a color-coded rectangle 
            color = motif_color_map.get(motif_seq, (0,0,0)) # setting tailored color according to motif_color_map, black if unrecognized / default 
            context.set_source_rgb(*color) # set MOTIF color = motif_color_map encoded
            context.rectangle(start, sequence_line_y - offset - motif_height - motif_margin, end-start, motif_height) # set MOTIF shape = rectangle & specifying its coordinates 
                #start & end = scaled position coordinates for exon (x left-side)
                #sequence_line_y - offset - motif_height = y top-side (same height as EXONS) 
                #end - start = width of rectangle 
                #motif_height = height of rectangle (same height as EXONS) 
            #context.rectangle(start, 110 + offset, end-start, 10)
                #start & end = scaled position coordinates for exon (x left-side)
                #50 = y top-side (same height as EXONS) 
                #end - start = width of rectangle 
                #50 = height of rectangle (same height as EXONS) 
            context.fill() #fill rectangle w green color 

        # 7. Add a color-coded key for motif types. 
        # setting key/legend position
        key_x = self.width - 750
        key_y = sequence_line_y - 800
        context.set_source_rgb(0,0,0) # text = black color 
        context.set_font_size(50) # font size 
        # drawing key 
        for motif, color in motif_color_map.items(): # for every pair in motif_color_map {}
            # draw a colored-box for each motif seq. 
            context.set_source_rgb(*color) # tailor color for each pair (motif_seq)
            context.rectangle(key_x,key_y,20,50) # draw a small colored box / pair 
            context.fill() 
            #draw motif name next to each colored box 
            context.set_source_rgb(0,0,0) # Black text 
            context.move_to(key_x+30,key_y+35) # position labels next to each color coded box 
            context.show_text(motif) # display motif seq string in legend 
            #move down for next entry in legend 
            key_y +=60 
        
        # 8. Add title 
        # 7. Add a color-coded key for motif types. 
        # setting key/legend position
        title_x = self.width/2 - 1000 # position @ ~ half way point 
        title_y = self.height/10 # position near top 
        context.move_to(title_x,title_y) # position title @ given coordinates
        context.set_source_rgb(0,0,0) # text = black color 
        context.set_font_size(75) # font size 
        # writing title 
        context.show_text(seq_title) # display sequence header as title

        # Save output graphic 
        return surface

# CLASS: COMBINE ALL VISUALIZATIONS --> one master .png file 
class Visualization_Combiner:
    def __init__(self): 
        self.visualizer = Visualizer() # allows use of the previous Visualizer() class 
        self.cairo_graphics = [] # init. empty list to store Visualizer() outputs (PyCairo Surfaces)

    def add_visualization(self, scaled_regions, scaled_motif_positions,seq_title): 
        "Creating Pycairo surface (via Visualizer.visualizer_seq()) per scaled sequence/motif coordinates & adding to a list (cairo_graphics[])"
        # Generate surface for given coordinates 
        surface = self.visualizer.visualize_seq(scaled_regions,scaled_motif_positions,seq_title)
        self.cairo_graphics.append(surface) # add surface to list cairo_graphics[]

    def combine_graphics(self,output_file): #output_file = name of output .png you want 
        "Combines all graphics from cairo_graphics [] into one master .png file"
        # Set dimensions for a larger surface to hold all combined graphics / seq. visualizations 
        total_width = max(surface.get_width() for surface in self.cairo_graphics) # total width = maximum width out of all cairo_graphics' widths 
        total_height = sum(surface.get_height() for surface in self.cairo_graphics) # total height = sum of all cairo_graphics' heights 

        combined_surface = cairo.ImageSurface(cairo.FORMAT_RGB24, total_width, total_height) # setting combined graphic dimensions & RGB range (the canvas) 
        combined_context = cairo.Context(combined_surface) # setting up the combined "paintbrush" for PyCairo 

        # Actually put each grphic onto larger surface (.png) space 
        current_y = 0 #init. empty variable to track where to paste each new graphic encountered in the loop (starting w/ leftmost edge)
        #loop through each graphic in cairo_graphics [] 
        for surface in self.cairo_graphics: 
            combined_context.set_source_surface(surface,0,current_y) # tell PyCairo to draw current surface onto the combined surface (combined_context) @ the bottom (0) @ the starting position (current_y) 
            combined_context.paint() # ACTUALLY paste the current surface onto previously specific position 
            current_y+= surface.get_height() # increment starting position for next surface (graphic) according to the width of the current surface  
        # Save the finalized combined surface/graphic onto a .png! 
        combined_surface.write_to_png(output_file)



# MAIN FUNCTION: connecting all the classes together to carry out their functions. 
def main(): 
# 1. Create Sequence Objects 
    # Create a Sequence object
    sequence_obj = Sequence('','') # setting (header, sequence) as empty strings ('') for now. The object will know to parse these out from the FASTA file when run below.
    # FASTA -- [parse_fasta]---> list of seq. obj.s
    sequences = sequence_obj.parse_fasta(input_fasta_file) 
# 2. Read Motifs from File 
    # MOTIF.txt --> motifs (list of motifs) 
    with open(input_motifs_file, 'r') as file: 
        motifs = file.readlines() # motifs = list containing each motif (line) listed in the input motif file 
    # Create Motif objects --> store in motif_objects [] 
    motif_objects = [Motif(motif.strip()) for motif in motifs] # loop throgh motifs [] --> strips new line chr. --> creates an obj. for each one 
# 3. Find (A) Motifs & (B) Exons & Introns in each given Sequence 
    combined_visualizer = Visualization_Combiner() # init. combined_visualizer
    seq_tracker=0
    for sequence in sequences: # loop through sequence objects 
        # Init. lists to hold all motif positions & all exon_intrion positions
        all_motif_positions = [] 
        all_exon_intron_positions = [] 
    #(A) Search for Motifs in Seq. --> Store coordinates in all_motif_positions []
        for motif_obj in motif_objects: # loop through motif objects PER sequence object 
            motif_positions = motif_obj.find_motif_in_seq(sequence.sequence) # search for motif in sequence (specified by sequence(obj. name).sequence(actual seq. string stored in obj.)) ---> coordinates 
            all_motif_positions.extend(motif_positions) # append the coordinates to the init.ed list each time 
    #(B) Search for Exon & Introns in Seq. --> Store coordinates in all_exon_intron_positions[]
        exon_intron_regions = motif_obj.find_exons_introns_in_seq(sequence.sequence) # search for exons/introns in sequence (specified by sequence(obj. name).sequence(actual seq. string stored in obj.)) ---> coordinates 
        all_exon_intron_positions.extend(exon_intron_regions) # append the coordinates to the init.ed list each time 
    # 4. Visualize the Sequence & its Features (exons, introns, motifs)
        # Create Visualizer object to handle the seq. visualization
        visualizer_obj = Visualizer() # (inputs were already predefined for class)
        # Scaling coordinates to fit dimensions of graphic: all_exon_intron_positions [] + all_motif_positions [] --[scaling_feature_coordinates]--> scaled_regions + scaled_motifs
        scaled_regions, scaled_motifs = visualizer_obj.scaling_feature_coordinates(all_exon_intron_positions,all_motif_positions)
        # Draw re-scaled coordinates via [visualize_seq]
        seq_title=sequence.header
        combined_visualizer.add_visualization(scaled_regions,scaled_motifs,seq_title) #running visualizer_obj.visualize_seq for every single seq. & adding graphic to a list 
    # 5. Combine all visualizations onto one master .png 
    combined_visualizer.combine_graphics(output_png_file)
    print(f"Master Visualization saved as {output_png_file}")

    # DEBUGGING NOTEBOOK 
    #print(sequences[0])
    #print(sequences[1])
    #print(sequences[2])
    #print(sequences[3])
    #print(all_motif_positions)
    #print(all_exon_intron_positions)
    #print(scaled_motifs)
    #print(scaled_regions)


# RUNNING MAIN FUNCTION: 
    # writing this as an ifelse ensures that main() fxn is only called if this script is directly executed (NOT when imported as a module to be utilized in other contexts)
if __name__=='__main__':
   main()