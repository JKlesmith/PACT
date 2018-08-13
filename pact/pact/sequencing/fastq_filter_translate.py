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

"""fastq_filter_translate - Filter and translate merged FASTQ files"""

from sys import version_info, maxsize

#Check to see if we're on 64-bit
if maxsize < 2**32:
    print("[PACT WARNING] This is not a 64-bit system, memory issues may occur with large files!")

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[FASTQ Filter/Translate Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from collections import Counter
from itertools import islice
from math import pow
from multiprocessing import Pool as ThreadPool
from multiprocessing import cpu_count
from statistics import median, mean
from pact.pact_common import file_checker

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

#Multiprocessing doesn't work inside of classes (workaround to expose the method as a top level function)
def unwrap_self_translate_dna_anchors(arg, **kwarg):
    return fastq_filter_translate.filter_dna_anchors(*arg, **kwarg)

def unwrap_self_translate_dna_nonanchored(arg, **kwarg):
    return fastq_filter_translate.filter_dna_nonanchored(*arg, **kwarg)

def unwrap_self_dna_seq_align_worker(arg, **kwarg):
    return fastq_filter_translate.dna_seq_align_worker(*arg, **kwarg)

class fastq_filter_translate:
    """Read and process FASTQ files"""

    def __init__(self, settings_dict):
        """Setup our class"""

        #Assign our varibles to the class
        self.wtdna = settings_dict['WTDNA']
        self.wtaa = settings_dict['WTAA']
        self.FirstAAMutated = int(settings_dict['FirstAAMutated'])-1
        self.LastAAMutated = int(settings_dict['LastAAMutated'])-1

        self.WTAARegion = self.wtaa[self.FirstAAMutated:self.LastAAMutated + 1]
        self.WTDNARegion = self.wtdna[self.FirstAAMutated*3:(self.LastAAMutated + 1) * 3]

        self.WTDNALen = len(self.WTDNARegion)
        self.WTAALen = len(self.WTAARegion)

        self.fiveprime = settings_dict['5pAnchor']

        #Set the q score limits
        try:
            self.q_average = pow(10, (-1 * int(settings_dict['QAverage'])) / 10)
            self.q_limit = pow(10, (-1 * int(settings_dict['QLimit'])) / 10)
        except ValueError:
            self.q_average = pow(10, (-20 / 10))
            self.q_limit = pow(10, (-10 / 10))

        #Check to see if the number of processes is logical
        if int(settings_dict['Processes']) <= 0 or settings_dict['Processes'] is None:
            self.processes = cpu_count()
        else:
            self.processes = int(settings_dict['Processes'])

        #Check to see if the mut count threshold is logical
        if int(settings_dict['MutThreshold']) < 1:
            self.mutcountthreshold = 1
        else:
            self.mutcountthreshold = int(settings_dict['MutThreshold'])

        #Check to see if we will enable anchors
        if settings_dict['Enable_Anchors'].lower() == "false":
            self.enable_anchors = False
        else:
            self.enable_anchors = True

        #Check to see if the fastq file exists
        if file_checker(settings_dict['In_File']):
            self.in_file = settings_dict['In_File']
        
        #Get the output prefix
        self.out_prefix = settings_dict['Out_Prefix']

        #Setup the translation table
        self.translation_table = {'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
                        'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
                        'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
                        'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
                        'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
                        'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
                        'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
                        'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
                        'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
                        'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
                        'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
                        'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
                        'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
                        'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
                        'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
                        'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'}

        return

    def dna_seq_align_worker(self, sequence):
        """Perform a non-gapping pairwise alignment."""

        #sequence_list = [Smaller string that we are searching, Larger string we are matching to]
        #If the anchors are not usable on the reads let's predict what the anchors would be through alignments (very computationally intensive)
        #The FASTQ is the columns as we are matching to it (as we are matching the translate WT region)

        #Get the length of the WT and search
        len_fastq = len(sequence)

        #Create a matrix with value 0 with WT length as the number of cols and the search as the number of rows
        list_matrix = [[1 if self.WTDNARegion[y] == sequence[x] else 0 for x in range(len_fastq)] for y in range(self.WTDNALen)]

        #Skew the matrix such that the diagnols are now columns
        for i in range(0, self.WTDNALen): #Rows = WT Letters
            list_matrix[i][0:0] = [0] * (self.WTDNALen - i - 1) #Add 0 (len(search)-1) to the front of the row
            list_matrix[i][len(list_matrix[i]):len(list_matrix[i])] = [0] * (i) #Add 0 i times to the end of the row

        #Get the column sums
        max_value = 0
        max_index = 0
        for j in range(0, len_fastq + self.WTDNALen - 1): #Columns
            sum_column = 0
            for i in range(0, self.WTDNALen): #Rows
                sum_column = sum_column + list_matrix[i][j]
            if sum_column > max_value:
                max_index = j
                max_value = sum_column

        #Return the predicted 5' anchor
        return sequence[max_index - self.WTDNALen + 1 - 9:max_index - self.WTDNALen + 1]

    def anchor_predict(self, sequence, mode_predict=True):
        """Predict the anchors using a non-gapping msa and then verify."""

        #Cap the number of reads to use
        if len(sequence) >= 1000:
            reads_to_sample = 1000
            reads_step = 10
        elif len(sequence) >= 100 and len(sequence) < 1000:
            reads_to_sample = len(sequence)
            reads_step = 5
        else:
            reads_to_sample = len(sequence)
            reads_step = 2

        #Setup a string that we can add to
        string_ret = ""

        #If we are predicting our anchors
        if mode_predict:
            print("[FASTQ Filter/Translate] Predicting anchor from multiple string alignment method.")

            #Make a list of sequences to test
            seq_to_align = [sequence[i][0] for i in range(0, reads_to_sample, reads_step) if 'N' not in sequence[i][0]]

            #Make a pool of worker threads
            pool = ThreadPool(self.processes)
            results = pool.map(unwrap_self_dna_seq_align_worker, zip([self]*len(seq_to_align), seq_to_align))
            pool.close()
            pool.join()

            #Let's get the most common sequence then assign the anchors to the class
            counts5p = Counter(result for result in results)
            self.fiveprime = counts5p.most_common(1)[0][0]

            string_ret = ("[FASTQ Filter/Translate] Anchor information:\n"
            "[FASTQ Filter/Translate] Sequences aligned: " + str(len(seq_to_align)) + "\n"
            "[FASTQ Filter/Translate] Predicted 5' Anchor: " + counts5p.most_common(1)[0][0] + "\tCounts: " + str(counts5p.most_common(1)[0][1]))

        #Final sanity check
        if len(self.fiveprime) == 0:
            print("[FASTQ Filter/Translate Error] Anchor length error.")
            quit()

        #Check the average read start of our 5' anchor
        self.median_start = int(median([sequence[i][0].find(self.fiveprime) for i in range(0, reads_to_sample, reads_step)]))
        
        #Length of anchor
        self.len_anchor = len(self.fiveprime)

        return string_ret

    def filter_dna_anchors(self, sequence):
        """Filter the FASTQ read (Called by FASTQFilter)"""

        #Calculate the length of the sequence using the anchor
        read_start = sequence[0].find(self.fiveprime)
        
        #If we can not find our anchor (could be a point mutation in the 5')
        if read_start == -1:
            tlsequence = sequence[0][self.median_start + self.len_anchor:self.median_start + self.len_anchor + self.WTDNALen]
        else:
            tlsequence = sequence[0][read_start + self.len_anchor:read_start + self.len_anchor + self.WTDNALen]

        #If it can not find the anchor (i.e. it has a mutation then let's reject it)
        #Using the anchors let's see if it matches the wild-type
        if len(tlsequence) != self.WTDNALen:
            return "length"

        #Reject sequences with N base calls
        if 'N' in tlsequence:
            return "ncall"

        #Convert the string to ascii then to the Q score probability using the Illumina 1.8 encoding
        qual_sequence = sequence[1][read_start + self.len_anchor:read_start + self.len_anchor + self.WTDNALen]
        scorelist = list(pow(10, (-1 * (int(code)-33)) / 10) for code in bytearray(qual_sequence, 'ascii'))

        if mean(scorelist) > self.q_average or max(scorelist) > self.q_limit:
            return "quality"

        #Let's translate our sequence
        seq_aa = ''.join([self.translation_table[tlsequence[i:i+3]] for i in range(0, self.WTDNALen, 3)])

        #Count the number of mismatches to the wild-type
        count = sum(1 for a, b in zip(seq_aa, self.WTAARegion) if a != b)

        #Wild-Type
        if count == 0:
            return ["WT\tWT", seq_aa, tlsequence, count]

        #Calculate where and what mutations are in the sequence
        seqids = ""
        muts = ""
        for i in range(0, self.WTAALen):
            if seq_aa[i] != self.WTAARegion[i]:
                seqids = ''.join([seqids, str(i), ","])
                muts = ''.join([muts, seq_aa[i], ","])

        #Return non-wildtype mutations
        return [''.join([seqids.rstrip("\n\r,"), "\t", muts.rstrip(",")]), seq_aa, tlsequence, count]

    def filter_dna_nonanchored(self, sequence):
        """Filter the FASTQ read (Called by FASTQFilter)"""
       
        #Calculate the length of the sequence
        tlsequence = sequence[0]
        tlread_len = len(sequence[0])

        #Reject if our read length does not match the wild-type
        if tlread_len != self.WTDNALen:
            return "length"

        #Reject sequences with N base calls
        if 'N' in tlsequence:
            return "ncall"

        #Calculate the quality scores
        scorelist = list(pow(10, (-1 * (int(code)-33)) / 10) for code in bytearray(sequence[1], 'ascii'))

        if mean(scorelist) > self.q_average or max(scorelist) > self.q_limit:
            return "quality"

        #Let's translate our sequence
        seq_aa = ''.join([self.translation_table[tlsequence[i:i+3]] for i in range(0, tlread_len, 3)])

        #Count the number of mismatches to the wild-type
        count = sum(1 for a, b in zip(seq_aa, self.WTAARegion) if a != b)

        #Wild-Type
        if count == 0:
            return ["WT\tWT", seq_aa, tlsequence, count]

        #Calculate where and what mutations are in the sequence
        seqids = ""
        muts = ""
        for i in range(0, self.WTAALen):
            if seq_aa[i] != self.WTAARegion[i]:
                seqids = ''.join([seqids, str(i), ","])
                muts = ''.join([muts, seq_aa[i], ","])

        #Return non-wildtype mutations
        return [''.join([seqids.rstrip("\n\r,"), "\t", muts.rstrip(",")]), seq_aa, tlsequence, count]

    def fastq_filter_translate(self):
        """Parse the FASTQ file"""

        #Update the user
        print("[FASTQ Filter/Translate] Reading in the FASTQ file: "+self.in_file)

        #Read the FASTQ file into a list
        with open(self.in_file, 'r') as input_fastq:
            dnaseq = list(islice(input_fastq, 1, None, 2))
        
        #Combine the DNA sequence and quality into one
        list_read_quality = [[dnaseq[i].rstrip('\n\r'), dnaseq[i+1].rstrip('\n\r')] for i in range(0, len(dnaseq), 2)]

        #Check for the correct anchors
        str_anchorret = ""
        if self.enable_anchors:
            if (self.fiveprime == "" or self.fiveprime == '\"\"' or len(self.fiveprime) == 0):
                print("[FASTQ Filter/Translate] DNA anchor was not provided, running prediction.")
                str_anchorret = self.anchor_predict(list_read_quality)
            else:
                print("[FASTQ Filter/Translate] Testing provided anchor.")
                str_anchorret = self.anchor_predict(list_read_quality, False)

        #Update the user
        print("[FASTQ Filter/Translate] Mapping processes to filter and translate the DNA.")

        #Make a pool of worker threads
        pool = ThreadPool(self.processes)

        #Open the DNA sequences in their own threads and return the results
        try:
            if self.enable_anchors:
                results = pool.map(unwrap_self_translate_dna_anchors, zip([self]*len(list_read_quality), list_read_quality))
            else:
                results = pool.map(unwrap_self_translate_dna_nonanchored, zip([self]*len(list_read_quality), list_read_quality))
        except MemoryError:
            print("[FASTQ Filter/Translate Error] Out of memory. Please use a 64bit install of python or install more ram.")
            quit()

        #Close the pool and make sure the results are in
        pool.close()
        pool.join()

        #Update the user
        print("[FASTQ Filter/Translate] Reducing the results and writing the output.")

        #Counters for the FASTQRead Functions
        count_fail_length = 0
        count_fail_nbase = 0
        count_fail_quality = 0
        count_wildtype = 0
        count_above_threshold = 0
        count_within_threshold = 0
        count_total_recorded = 0

        #Open the FASTQ read file
        output_fastq = open(self.out_prefix+'_Read.tsv', 'w')

        #Add the header to the file
        output_fastq.write("Location\tMutations\tAA_Sequence\tDNA_Sequence\tNum_Muts\n")

        #Go through each result
        for result in results:
            
            #Check to see if the length
            if result == "length":
                count_fail_length = count_fail_length + 1
                continue

            #Check to see if the nbase filter passed
            if result == "ncall":
                count_fail_nbase = count_fail_nbase + 1
                continue

            #Check to see if the quality filter passed
            if result == "quality":
                count_fail_quality = count_fail_quality + 1
                continue

            #Check to see which type of mutation it is
            if result[3] == 0:
                count_wildtype = count_wildtype + 1
            elif result[3] > 0 and result[3] <= self.mutcountthreshold:
                count_within_threshold = count_within_threshold + 1
            else:
                count_above_threshold = count_above_threshold + 1

            #Count the total number of results
            count_total_recorded = count_total_recorded + 1

            #Write to our file
            output_fastq.write(''.join([result[0], "\t", result[1], "\t", result[2], "\t", str(result[3]), "\n"]))

        #Close the file
        output_fastq.close()

        #Print information
        output_string = str_anchorret + "\n" + '\n'.join([
        "[FASTQ Filter/Translate] FASTQ Filename: " + self.in_file,
        "[FASTQ Filter/Translate] Sequences rejected due to incorrect DNA length relative to WT: " + str(count_fail_length),
        "[FASTQ Filter/Translate] Sequences rejected due to having N called bases: " + str(count_fail_nbase),
        "[FASTQ Filter/Translate] Sequences rejected due to poor read quality (Q average or Q lower limit): " + str(count_fail_quality),
        "[FASTQ Filter/Translate] Sequences passed QC and recorded: " + str(count_total_recorded),
        "[FASTQ Filter/Translate] Sequences with 0 nonsynonymous mutations: " + str(count_wildtype),
        "[FASTQ Filter/Translate] Sequences with 1 to " + str(self.mutcountthreshold) + " nonsynonymous mutations: " + str(count_within_threshold),
        "[FASTQ Filter/Translate] Sequences with more than " + str(self.mutcountthreshold) + " nonsynonymous mutations: " + str(count_above_threshold),
        ])

        #Update the user
        print(output_string)

        return output_string

if __name__ == '__main__':

    #Remind the user that the module needs to be ran within the context of PACT
    print("[FASTQ Filter/Translate Error] This module needs to be ran within the context of PACT.")
