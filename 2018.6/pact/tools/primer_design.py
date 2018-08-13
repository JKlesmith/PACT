#!/usr/bin/env python3
#
# Protein Analysis and Classifier Toolkit
# Copyright (C) 2018 Justin R. Klesmith
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""primer_design - design SSM oligos for Nicking SSM"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Primer Design Error] Your Python interpreter is too old. Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from math import exp, log
from multiprocessing import Pool as ThreadPool

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.5"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

#Multiprocessing doesn't work inside of classes (workaround to expose the method as a top level function)
def unwrap_self_make_oligo(arg, **kwarg):
    return primer_design.make_oligo(*arg, **kwarg)

class primer_design:
    """Design primers for nicking SSM"""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""

        #Assign our varibles to the class
        self.wtdna = settings_dict['WTDNA'].upper()

        #Check to see if the number of processes is logical
        if int(settings_dict['Processes']) <= 0 or settings_dict['Processes'] is None:
            from multiprocessing import cpu_count
            self.processes = cpu_count()
        else:
            self.processes = int(settings_dict['Processes'])

        #Import the mutation design list
        self.list_mutation_design = []
        for group in literal_eval(settings_dict['mutcodons']):
            #If we have n between two numbers
            if 'n' in group:
                #We need exactly three to expand the range
                if len(group) == 3:
                    #Get the start and end position, then iterate through the two points
                    list_expandgroup = []
                    for i in range(group[0], group[-1] + 1):
                        list_expandgroup.append(i)
                    self.list_mutation_design.append(list_expandgroup)
                else:
                    print("[Error] The codon list is missing either a start or end point around the 'n' location.")
                    exit(1)
            else:
                #Add the other ranges without changing them
                self.list_mutation_design.append(group)

        #Set the output file prefix
        self.outfile_prefix = settings_dict['Out_Prefix']

        return

    def sequence_params(self, sequence, mismatch = 3):
        """Calculate the length, gc%, mismatch% and Tm of the sequence."""

        #Make sure it is in upper case
        sequence = sequence.upper().rstrip('\n')

        #Calculate the length of the sequence
        sequence_length = len(sequence)

        #Return the percent of GC in the sequence as a whole number
        sequence_gc = int(((sequence.count("C") + sequence.count("G")) / sequence_length) * 100)

        #Return the mismatch percent as a whole number
        sequence_mismatch = int((mismatch / sequence_length) * 100)

        #Calculate the Tm
        sequence_tm = round(81.5 + 0.41*sequence_gc - (675/sequence_length) - sequence_mismatch, 2)

        return sequence_length, sequence_gc, sequence_tm

    def sequence_params_sides(self, sequence_gapped):
        """Return the side sequences and the length."""

        #Then calculate the gapped sequence (with NNK)
        #Get the location of the degenerate codon
        nnk_index = sequence_gapped.upper().find("NNK")
        nnc_index = sequence_gapped.upper().find("NNC")
        nnn_index = sequence_gapped.upper().find("NNN")

        #Only use the index if it is there
        if nnk_index != -1:
            deg_codon_index = nnk_index
        elif nnn_index != -1:
            deg_codon_index = nnn_index
        elif nnc_index != -1:
            deg_codon_index = nnc_index #Set as last as NNNC will fall in here and mess things up

        #Return the sequence to the left and right of the codon
        seq_5p = sequence_gapped[:deg_codon_index].upper().rstrip('\n')
        seq_3p = sequence_gapped[deg_codon_index + 3:].upper().rstrip('\n')
        
        #Get the sequence length
        seq_5p_len = len(seq_5p)
        seq_3p_len = len(seq_3p)
        
        return seq_5p, seq_3p, seq_5p_len, seq_3p_len

    def nn_energy(self, sequence_nnk, sequence_wt):
        """Calculate the free energy and cost using nearest neighbor model (SantaLucia)"""

        #Clean up our inputs
        sequence_nnk = sequence_nnk.upper().rstrip('\n')
        sequence_wt = sequence_wt.upper().rstrip('\n')

        #Setup the lookup tables for dS and dH
        ds = {
            "AA": -22.2,
            "TT": -22.2,
            "AT": -20.4,
            "TA": -21.3,
            "CA": -22.7,
            "TG": -22.7,
            "GT": -22.4,
            "AC": -22.4,
            "CT": -21,
            "AG": -21,
            "GA": -22.2,
            "TC": -22.2,
            "CG": -27.2,
            "GC": -24.4,
            "GG": -19.9,
            "CC": -19.9
            }

        ds_init_termbase = {
            "G":-2.8,
            "C":-2.8,
            "A":4.1,
            "T":4.1
            }

        dh = {
            "AA": -7.9,
            "TT": -7.9,
            "AT": -7.2,
            "TA": -7.2,
            "CA": -8.5,
            "TG": -8.5,
            "GT": -8.4,
            "AC": -8.4,
            "CT": -7.8,
            "AG": -7.8,
            "GA": -8.2,
            "TC": -8.2,
            "CG": -10.6,
            "GC": -9.8,
            "GG": -8,
            "CC": -8
            }

        dh_init_termbase = {
            "G":0.1,
            "C":0.1,
            "A":2.3,
            "T":2.3
            }

        #for dg 68oC = 341.1K, 37oC = 310.2K dG = dH - TdS
        #First calculate the free energy of the entire ungapped template
        list_dh_wt = list(dh[sequence_wt[i:i+2]] for i in range(0, len(sequence_wt)) if len(sequence_wt[i:i+2]) == 2)
        list_ds_wt = list(ds[sequence_wt[i:i+2]] for i in range(0, len(sequence_wt)) if len(sequence_wt[i:i+2]) == 2)
        
        wt_h = sum(list_dh_wt) + dh_init_termbase[sequence_wt[-1:]]
        wt_s = sum(list_ds_wt) + ds_init_termbase[sequence_wt[-1:]]

        wt_g = round(wt_h - (341.1 * wt_s * 0.001), 2) #0.001 as it's in cal and needs to be converted to kcal

        #Calculate the sides
        seq_5p, seq_3p, seq_5p_len, seq_3p_len = self.sequence_params_sides(sequence_nnk)

        #Calculate the dH and dS
        list_dh_5p = list(dh[seq_5p[i:i+2]] for i in range(0, seq_5p_len) if len(seq_5p[i:i+2]) == 2)
        list_ds_5p = list(ds[seq_5p[i:i+2]] for i in range(0, seq_5p_len) if len(seq_5p[i:i+2]) == 2)

        list_dh_3p = list(dh[seq_3p[i:i+2]] for i in range(0, seq_3p_len) if len(seq_3p[i:i+2]) == 2)
        list_ds_3p = list(ds[seq_3p[i:i+2]] for i in range(0, seq_3p_len) if len(seq_3p[i:i+2]) == 2)

        gap_h = sum(list_dh_5p) + dh_init_termbase[seq_5p[-1:]] + sum(list_dh_3p) + dh_init_termbase[seq_3p[-1:]]
        gap_s = sum(list_ds_5p) + ds_init_termbase[seq_5p[-1:]] + sum(list_ds_3p) + ds_init_termbase[seq_3p[-1:]]
        
        gap_g = round(gap_h - (341.1 * gap_s * 0.001), 2)

        #Energy cost of mismatches is (Eungapped - Egapped)/Eungapped * 100
        cost = round(((wt_g - gap_g)/wt_g)*100, 2)

        #Return the dG ungapped, dG gapped, and energy cost %
        return wt_g, gap_g, cost

    def phusion_tm(self, sequence_nnk):
        """Calculate the Tm in a phusion reaction
        Sources: 
        NEB website to which equation to use (Breslauer)
        SigmaAldrich, MuPlex on the helix initiation correction -10.8 cal/Mol (*2 for the NNK oligo, treated as two oligos)
        Buffer cation concentration is from a testset against the NEB website tool
        Adjusted for Nicking primer concentration of 0.0003524uM
        """

        #Clean up our inputs
        sequence_nnk = sequence_nnk.upper().rstrip('\n')

        #Setup the lookup tables for dS and dH
        ds = {
            "AA": -24,
            "TT": -24,
            "AT": -23.9,
            "TA": -16.9,
            "CA": -12.9,
            "TG": -12.9,
            "GT": -17.3,
            "AC": -17.3,
            "CT": -20.8,
            "AG": -20.8,
            "GA": -13.5,
            "TC": -13.5,
            "CG": -27.8,
            "GC": -26.7,
            "GG": -26.6,
            "CC": -26.6
            }

        dh = {
            "AA": -9.1,
            "TT": -9.1,
            "AT": -8.6,
            "TA": -6,
            "CA": -5.8,
            "TG": -5.8,
            "GT": -6.5,
            "AC": -6.5,
            "CT": -7.8,
            "AG": -7.8,
            "GA": -5.6,
            "TC": -5.6,
            "CG": -11.9,
            "GC": -11.1,
            "GG": -11,
            "CC": -11
            }

        #Calculate the sides
        seq_5p, seq_3p, seq_5p_len, seq_3p_len = self.sequence_params_sides(sequence_nnk)

        list_dh_5p = sum(list(dh[seq_5p[i:i+2]] for i in range(0, seq_5p_len) if len(seq_5p[i:i+2]) == 2))
        list_ds_5p = sum(list(ds[seq_5p[i:i+2]] for i in range(0, seq_5p_len) if len(seq_5p[i:i+2]) == 2))

        list_dh_3p = sum(list(dh[seq_3p[i:i+2]] for i in range(0, seq_3p_len) if len(seq_3p[i:i+2]) == 2))
        list_ds_3p = sum(list(ds[seq_3p[i:i+2]] for i in range(0, seq_3p_len) if len(seq_3p[i:i+2]) == 2))

        gap_h = list_dh_5p + list_dh_3p
        gap_s = list_ds_5p + list_ds_3p
        
        gap_tm = round(((gap_h * 1000)/(-21.6 + gap_s + (1.987 * log(3.524*10**-10/4)))) - 273.15 + 16.6 * log(0.28), 2)

        #Return the tm
        return gap_tm

    def qc_oligo_filter(self, sequence):
        """Do a first pass quality check of the oligo and reject with filters (expects oligo with NNK)."""

        #Make sure our sequence is in upper case and without the newline
        sequence = sequence.upper().rstrip('\n')

        #Hard limits from emperical data
        length_max = 60 #59 observed, 60 theoretical
        length_min = 23 #observed
        min_gc_percent = 30 #30 observed, 40 recommended
        max_gc_percent = 75 #observed
        min_melt_temp = 63 #observed
        max_melt_temp = 82 #observed

        #Return the oligo properties
        sequence_length, sequence_gc, sequence_tm = self.sequence_params(sequence, 3)

        #Check the min and max length
        if sequence_length > length_max:
            return False, "LENHIGH"
        elif sequence_length < length_min:
            return False, "LENLOW"

        #Check the min and max melt temp
        if sequence_tm > max_melt_temp:
            return False, "TMHIGH"
        elif sequence_tm < min_melt_temp:
            return False, "TMLOW"

        #Check the min and max gc content
        if sequence_gc > max_gc_percent:
            return False, "GCHIGH"
        elif sequence_gc < min_gc_percent:
            return False, "GCLOW"

        #Check to see if the GC is on the GC vs Len curve (fits a one phase exp decay curve at the edges max/min)
        #Measurements originated from oligos designed with Agilent QC program and successfully used with Nicking
        #Input the len of the sequence and output a expected GC%
        #The fits were performed in GraphPad Prism

        #Lower edge
        if sequence_gc < (185.7 * exp(-0.06797 * sequence_length)) + 21.65:
            return False, "GCFITLOW"

        #Upper edge
        if sequence_gc > (205.1 * exp(-0.05125 * sequence_length)) + 28.64:
            return False, "GCFITHIGH"

        #Get the information of the sequence on the sides of the NNK
        seq_5p, seq_3p, seq_5p_len, seq_3p_len = self.sequence_params_sides(sequence)
        p5_length, p5_gc, p5_tm = self.sequence_params(seq_5p, 0)
        p3_length, p3_gc, p3_tm = self.sequence_params(seq_3p, 0)

        #Check for the minimum length (min 6 observed but is rare)
        if p5_length < 10 or p3_length < 10:
            return False, "LENSIDELOW"

        #Check for the maximum length (max 37 observed but is rare, 39 is theoretical max if max len diff of sides is 18)
        if p5_length > 39 or p3_length > 39:
            return False, "LENSIDEHIGH"

        #Check for the minimum GC content (no need to check for max as obs was 100% at the max)
        if p5_gc < 21 or p3_gc < 21:
            return False, "GCSIDELOW"

        #Check for GC of the oligo on each side (this is related to above).
        #The GC vs Len can be weakly fit with a exp decay curve
        #Input the len and the expected GC is returned

        side_gc_5p_low = (148.9 * exp(-0.1929 * p5_length)) + 18.11
        side_gc_5p_high = (193.4 * exp(-0.03400 * p5_length)) - 24.26
        
        side_gc_3p_low = (148.9 * exp(-0.1929 * p3_length)) + 18.11
        side_gc_3p_high = (193.4 * exp(-0.03400 * p3_length)) - 24.26

        if p5_gc > side_gc_5p_high or p3_gc > side_gc_3p_high:
            return False, "GCSIDEFITHIGH"

        if p5_gc < side_gc_5p_low or p3_gc < side_gc_3p_low:
            return False, "GCSIDEFITLOW"

        #The sequence passes.
        return True, "PASS"

    def qc_oligo_score(self, list_sequence):
        """From a list of sequences pick the best one from the dataset from a sequence score."""

        #Pull out the gapped sequence our combined list of gap, wt
        sequence = list_sequence[0].upper().rstrip('\n')

        #Calculate the paramaters of the sequence
        sequence_length, sequence_gc, sequence_tm = self.sequence_params(sequence)
        seq_5p, seq_3p, seq_5p_len, seq_3p_len = self.sequence_params_sides(sequence)
        p5_length, p5_gc, p5_tm = self.sequence_params(seq_5p, 0)
        p3_length, p3_gc, p3_tm = self.sequence_params(seq_3p, 0)
        phusion_tm = self.phusion_tm(sequence)

        #Total Score
        score = 0

        #Check the full length oligo

        #Check if the last character is a G or C
        if sequence[-1:] == "G" or sequence[-1:] == "C":
            score = score + 0.5

        #Check to see if the full length is -1 SD and +1.5 SD from average (36 +/- 5) (longer oligos are often better)
        if sequence_length >= (36 - 5) and sequence_length <= (36 + 8):
            score = score + 1

        #Check to see if the full GC is -1 SD and +1.5 SD from average (51 +/- 7)
        if sequence_gc >= (51 - 7) and sequence_gc <= (51 + 10):
            score = score + 1

        #Check to see if the full tm is 1 SD from the mode (very left skewed) (75 +/- 3)
        if sequence_tm >= (78 - 3) and sequence_tm <= (78 + 3):
            score = score + 1

        #Check to see if the phusion tm is 1 SD from median (dead center of the population) (left skewed) (68 +/- 4)
        if phusion_tm >= (68 - 4) and phusion_tm <= (68 + 4):
            score = score + 1

        #Do a hard check for the calculated Tm in phusion buffer for Nicking (min 50 but set at 60 for safety)
        if phusion_tm < 60:
            score = -10


        #Now check the sides of the oligo

        #Check to see if the side length is 1 SD from average (16 +/- 4)
        if p5_length >= (16 - 4) and p5_length <= (16 + 4):
            score = score + 0.5

        #Check to see if the side GC is 1 SD from average (57 +/- 13)
        if p5_gc >= (57 - 13) and p5_gc <= (57 + 13):
            score = score + 0.5

        #Check to see if the side length is 1 SD from average (17 +/- 3)
        if p3_length >= (17 - 3) and p3_length <= (17 + 3):
            score = score + 0.5

        #Check to see if the side GC is 1 SD from average (56 +/- 12)
        if p3_gc >= (56 - 12) and p3_gc <= (56 + 12):
            score = score + 0.5

        #Check to see if the side len diff is 1 SD from 0, actual is (-1 +/- 4.6)
        if abs(p5_length - p3_length) >= (0 - 4) and abs(p5_length - p3_length) <= (0 + 4):
            score = score + 2

        return [score, sequence, list_sequence[1], sequence_length, sequence_gc, round(sequence_tm, 2), round(phusion_tm, 2)]

    def make_oligo(self, codon_number):
        """Make all possible oligos and then return the best one for a given position."""

        #For the sides of the NNK, the min and max sequence to use
        side_len_min = 10
        side_len_max = 39

        #Compute the actual max based on the codon_number
        len_wtdna = len(self.wtdna)
        len_p5_max = (codon_number * 3) - 3 #The amount of sequence 5' of the codon
        len_p3_max = len_wtdna - (codon_number * 3) #The amount of sequence 3' of the codon

        #Reassign the max based on the above values
        if len_p5_max > side_len_max:
            len_p5_max = side_len_max

        if len_p3_max > side_len_max:
            len_p3_max = side_len_max

        #Assuming that the DNA is in frame, get the max sequence to the left and right of the codon
        side_5p = self.wtdna[codon_number * 3 - len_p5_max - 3:codon_number * 3 - 3]
        side_3p = self.wtdna[codon_number * 3:codon_number * 3 + len_p3_max]

        #Find the wild-type codon
        wild_type_codon = self.wtdna[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

        #Check to see if we have a sequence or we ran onto the plasmid
        #Will allow if the length is > 2 SD of average side length (24 bases)
        if len(side_5p) == 0 or len(side_5p) != 39:
            if len(side_5p) > 24:
                print("[Error] Missing 5' sequence. Will run as the length is >2 SD of avg side length. Codon #:"+str(codon_number))
                print("[Error] Have: " + str(len(side_5p)) + " bases, Need: 39 bases, Lower limit is 24 bases.")         
            else:
                print("[Error] Missing 5' sequence. Add 5' vector sequence and re-run. Codon #:"+str(codon_number))
                print("[Error] Have: " + str(len(side_5p)) + " bases, Need: 39")
                return [codon_number, None]

        if len(side_3p) == 0 or len(side_3p) != 39:
            if len(side_3p) > 24:
                print("[Error] Missing 3' sequence. Will run as the length is >2 SD of avg side length. Codon #:"+str(codon_number))
                print("[Error] Have: " + str(len(side_3p)) + " bases, Need: 39 bases, Lower limit is 24 bases.")         
            else:
                print("[Error] Missing 3' sequence. Add 3' vector sequence and re-run. Codon #:"+str(codon_number))
                print("[Error] Have: " + str(len(side_3p)) + " bases, Need: 39")
                return [codon_number, None]

        #Generate all of the 5' and 3' sequences
        list_5p = list(side_5p[-i:] for i in range(side_len_min, len_p5_max + 1))
        list_3p = list(side_3p[:i] for i in range(side_len_min, len_p3_max + 1))

        #Generate all of the oligo combinations
        list_oligos = list(list([seq5 + "NNK" + seq3, seq5 + wild_type_codon + seq3]) for seq5 in list_5p for seq3 in list_3p)
       
        #Filter each oligo
        list_filtered = []
        list_rejected = []
        for oligo in list_oligos:
            filter, status = self.qc_oligo_filter(oligo[0])

            #If the oligo passes the filter
            if filter:
                list_filtered.append(oligo)
            else:
                if status != "LENLOW" and status != "LENHIGH": #Most often it's a Tm or a GC issue
                    list_rejected.append(oligo)

        #Error out if there is no results (total filter rejection)
        if len(list_filtered) == 0:
            print("[Error] The hard filter has rejected all sequences. Codon: " + str(codon_number))
            print("The rejected sequences will be scored. However, it is suggested to check this oligo with a different tool.")
            
            #Score each oligo on our metrics
            list_scores = list(self.qc_oligo_score(oligo) for oligo in list_rejected)   

        else:
            #Score each oligo on our metrics
            list_scores = list(self.qc_oligo_score(oligo) for oligo in list_filtered)     

        #Sort the list for the top 15 scored sequences
        sorted_list_scores = sorted(list_scores, key=lambda x: -x[0])[:15]
        
        #Score each oligo on free energy and join the two lists
        list_energy = list(scored + list(self.nn_energy(scored[1], scored[2])) for scored in sorted_list_scores)

        #Sort the list for the lowest energy cost from column 9 (% cost)
        sorted_list_energy = sorted(list_energy, key=lambda x: x[9])

        #Return the codon number and the first primer
        return [codon_number] + sorted_list_energy[0]

    def primer_design(self):
        """Main entry point for the primer design application."""

        #Make a pool of worker threads
        pool = ThreadPool(self.processes)

        #Open the DNA sequences in their own threads and return the results
        try:
            results = pool.map(unwrap_self_make_oligo, zip([self]*len(self.list_mutation_design[0]), self.list_mutation_design[0]))
        except MemoryError:
            print("[Primer Design Error] Out of memory. Please use a 64bit install of python or install more ram.")
            exit(1)

        #Close the pool and make sure the results are in
        pool.close()
        pool.join()

        #Output the results to a file
        with open(self.outfile_prefix + '_OligoDesign.csv', 'w') as outfile:

            outfile.write(','.join([
                    "Number", #Number
                    "Name",
                    "Designed_Sequence", #Designed_Sequence
                    "WT_Sequence", #WT_Sequence
                    "Program_Score", #Program_Score
                    "Length", #Length
                    "GC%", #GC%
                    "Tm (oC)", #Tm (oC)
                    "Tm (oC) Phusion",
                    "dG WT (kcal/mole)", #dG Ungapped (kcal/mole)
                    "dG Design (kcal/mole)", #dG Gapped (kcal/mole)
                    "% Energy Cost" #% Energy Cost
                    ]) + "\n")

            #Loop the results
            for i in range(len(results)):

                #Skip mutations that we don't have
                if results[i][1] == None:
                    outfile.write(str(results[i][0]) + "," + str(results[i][0]) + "_NNK" + "\n")
                else:
                    outfile.write(','.join(map(str, [
                    results[i][0], #Number
                    str(results[i][0]) + "_NNK", #Number
                    results[i][2], #Designed Sequence
                    results[i][3], #WT Sequence
                    results[i][1], #Program Score
                    results[i][4], #Length
                    results[i][5], #GC%
                    results[i][6], #Tm (oC)
                    results[i][7], #Tm Phusion
                    results[i][8], #dG Ungapped (kcal/mole)
                    results[i][9], #dG Gapped (kcal/mole)
                    results[i][10], #% Energy Cost
                    ])) + "\n")

        return

if __name__ == '__main__':

    from sys import platform, modules
    from multiprocessing import freeze_support
    from argparse import ArgumentParser

    #If we are on windows we need to have freeze support for multithreading
    if platform.startswith('win'):
        #Not required if called from __main__ but the modules may be imported from a different script
        if __name__ != "__main__":
            modules["__main__"] = modules[__name__]
        freeze_support()

    #Parse the command line for inputs
    parser = ArgumentParser(description='Primer Design')
    parser.add_argument('-d', dest='dnasequence', action='store', help='DNA Sequence')
    parser.add_argument('-j', dest='processes', action='store', help='Number of processes to spawn')
    parser.add_argument('-o', dest='output_prefix', action='store', help='Output file name prefix')
    parser.add_argument('-m', dest='mutated_codons', action='store', help='What type of codon to mutate to')
    args = parser.parse_args()

    #Print a header
    print("[Primer Design] Primer Design\n"
            "[Primer Design] Author: "+__author__+"\n"
            "[Primer Design] Contact: "+__email__[0]+", "+__email__[1]+"\n"
            "[Primer Design] "+__copyright__+"\n"
            "[Primer Design] Version: "+__version__+"\n"
            "[Primer Design] License: "+__license__+"\n"
            "[Primer Design] Github [user: JKlesmith] (www.github.com)")

    pd_settings = {
        'WTDNA':args.dnasequence,
        'Processes':args.processes,
        'Out_Prefix':args.output_prefix,
        'mutcodons':args.mutated_codons
        }

    #Create our object
    classpd = primer_design(pd_settings)
    classpd.primer_design()
