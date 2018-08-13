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

"""fastq_merge - Merge FASTQ files"""

from sys import version_info, maxsize

#Check to see if we're on 64-bit
if maxsize < 2**32:
    print("[PACT WARNING] This is not a 64-bit system, memory issues may occur with large files!")

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[FASTQ Merge Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from itertools import islice
from multiprocessing import Pool as ThreadPool
from multiprocessing import cpu_count
from statistics import mode, stdev, StatisticsError

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

#Multiprocessing doesn't work inside of classes (workaround to expose the method as a top level function)
def unwrap_self_merge_pairwise(arg, **kwarg):
    return fastq_merge.merge_pairwise(*arg, **kwarg)

def unwrap_self_merge_pairwise_scout(arg, **kwarg):
    return fastq_merge.merge_pairwise_scout(*arg, **kwarg)

class fastq_merge:
    """Merge FASTQ Reads"""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""
        #Check to see if the number of processes is logical
        if settings_dict['Processes'] is None:
            self.processes = cpu_count()
        else:
            if int(settings_dict['Processes']) <= 0:
                self.processes = cpu_count()
            else:
                self.processes = int(settings_dict['Processes'])

        #Check our files
        if self.file_checker(settings_dict['Forward_FASTQ']):
            self.forward_fastq = settings_dict['Forward_FASTQ']

        if self.file_checker(settings_dict['Reverse_FASTQ']):
            self.reverse_fastq = settings_dict['Reverse_FASTQ']

        self.out_prefix = settings_dict['Out_Prefix']
        
        #Minimum coverage
        if settings_dict['Min_Coverage'] is None:
            self.min_coverage = 0.2
        else:
            if float(settings_dict['Min_Coverage']) < 0 or float(settings_dict['Min_Coverage']) > 1:
                self.min_coverage = 0.2
            else:
                self.min_coverage = float(settings_dict['Min_Coverage'])
        
        return

    def file_checker(self, file_location):
        """Check to see if the file exists"""
        from os.path import isfile

        try:
            if isfile(file_location):
                return True
            else:
                print("[File Import Error] The file: " + file_location + " cannot be found.")
                quit()
                return False
        except TypeError:
                print("[File Import Error] The file: " + str(file_location) + " cannot be found.")
                quit()
                return False

    def merge_pairwise_scout(self, line):
        """Return the index of the longest match."""

        #Define the dna and quality sequences
        #forward_dna = line[0][0]

        #Get the reverse comp of the 2nd read
        dict_rc = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        reverse_dna = ''.join([dict_rc[char] for char in line[1][0][::-1]])

        #Get the length of the forward and the reverse
        len_forward = len(line[0][0])

        #Set our counters and calculate the number of matches for wild-type
        max_value = sum(1 for a, b in list(zip(line[0][0], reverse_dna)) if a == b)
        max_index = 0
        float_length_fraction = max_value/len_forward
        direction = 1 #1 = FWD, 2 = REV

        #Iterate through each string combination
        for i in range(1, len_forward):

            #Calculate the number of matches in the forward read first orientation
            count = sum(1 for a, b in list(zip(line[0][0][-i:], reverse_dna[:i])) if a == b)
            
            #Assess if our count is greater than the running total
            if count > max_value:
                max_value = count
                max_index = i
                float_length_fraction = count/len_forward
                direction = 1

            #Calculate the number of matches in the reverse read first orientation
            count = sum(1 for a, b in list(zip(line[0][0][:i], reverse_dna[-i:])) if a == b)                

            #Assess if our count is greater than the running total
            if count > max_value:
                max_value = count
                max_index = i
                float_length_fraction = count/len_forward
                direction = 2

        return [max_index, float_length_fraction, direction]

    def merge_pairwise(self, line):
        """Merge the read pairs"""

        #Define the dna and quality sequences
        #forward_dna = line[0][0]
        #forward_quality = line[0][1]
        #reverse_quality = line[1][1][::-1]

        #Get the reverse comp of the 2nd read
        dict_rc = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        reverse_dna = ''.join([dict_rc[char] for char in line[1][0][::-1]])

        #Get the length of the forward and the reverse
        len_forward = len(line[0][0])

        #Do a length sanity check to make sure the read lengths are exactly the same
        if len_forward != len(reverse_dna):
            print("[FASTQ Merge Error] Read lengths are not the same.")
            quit()

        #Set our counters
        max_value = 0
        max_index = 0
        float_length_fraction = 0.0

        #First test our scouted list of indexes
        if self.direction_mode == 3:
            #If the initial scouting couldn't converage at a mode, then brute force it
            max_index, float_length_fraction, self.direction_mode = self.merge_pairwise_scout(line)
        else:
            #Else, use the guidance from the intial scouting
            for idx in self.list_tight_index:
                if idx == 0:
                    #Calculate the number of matches
                    count = sum(1 for a, b in list(zip(line[0][0], reverse_dna)) if a == b)
                else:
                    if self.direction_mode == 1:
                        #Calculate the number of matches
                        count = sum(1 for a, b in list(zip(line[0][0][-idx:], reverse_dna[:idx])) if a == b) 
                    elif self.direction_mode == 2:
                        count = sum(1 for a, b in list(zip(line[0][0][:idx], reverse_dna[-idx:])) if a == b)

                #Break our loop if we are at the max
                float_length_fraction = count/len_forward
                if float_length_fraction >= 0.7:
                    max_index = idx
                    break

        #If we fail then scout it out
        if float_length_fraction < 0.7:
            max_index, float_length_fraction, self.direction_mode = self.merge_pairwise_scout(line)

        #If our fraction is below the threshold match then reject
        if float_length_fraction < self.min_coverage:
            return [None]

        #Setup our DNA sequences and quality (based on the direction)
        if self.direction_mode == 1:
            #Our paired strings if the forward read is before the reverse
            forward_dna = ''.join([line[0][0], "-" * (len_forward - max_index)])
            reverse_dna = ''.join(["-" * (len_forward - max_index), reverse_dna])
            forward_quality = ''.join([line[0][1], "!" * (len_forward - max_index)])
            reverse_quality = ''.join(["!" * (len_forward - max_index), line[1][1][::-1]])
        else:
            #Our paired strings if the reverse is before the forward
            forward_dna = ''.join([reverse_dna, "-" * (len_forward - max_index)])
            reverse_dna = ''.join(["-" * (len_forward - max_index), line[0][0]])
            forward_quality = ''.join([line[1][1][::-1], "!" * (len_forward - max_index)])
            reverse_quality = ''.join(["!" * (len_forward - max_index), line[0][1]])

        #Forward Quality
        list_fwd_scores = list(int(code)-33 for code in bytearray(forward_quality, 'ascii'))

        #Reverse Quality
        list_rev_scores = list(int(code)-33 for code in bytearray(reverse_quality, 'ascii'))      

        #Loop the read paried list
        merged_read = ""
        merged_quality = ""
        for i in range(0, len(forward_dna)):

            #Check if there is a - in the reverse read (to build the start of the string)
            if forward_dna[i] != "-" and reverse_dna[i] == "-":
                merged_read = ''.join([merged_read, forward_dna[i]])
                merged_quality = ''.join([merged_quality, forward_quality[i]])
                continue

            #Check if the base call is the same
            if forward_dna[i] == reverse_dna[i]:
                #Add the base if they are the same (and take the higher q score)
                merged_read = ''.join([merged_read, forward_dna[i]])

                if list_fwd_scores[i] >= list_rev_scores[i]:
                    merged_quality = ''.join([merged_quality, forward_quality[i]])
                else:
                    merged_quality = ''.join([merged_quality, reverse_quality[i]])

                continue
            else:
                #Only take the higher of the two if the base call is not the same
                if list_fwd_scores[i] == list_rev_scores[i]:
                    merged_read = ''.join([merged_read, "N"])
                    merged_quality = ''.join([merged_quality, forward_quality[i]])
                elif list_fwd_scores[i] >= list_rev_scores[i]:
                    merged_read = ''.join([merged_read, forward_dna[i]])
                    merged_quality = ''.join([merged_quality, forward_quality[i]])
                else:
                    merged_read = ''.join([merged_read, reverse_dna[i]])
                    merged_quality = ''.join([merged_quality, reverse_quality[i]])

                continue

            #Check if there is a - in the forward read (to build the end of the string)
            if forward_dna[i] == "-" and reverse_dna[i] != "-":
                merged_read = ''.join([merged_read, reverse_dna[i]])
                merged_quality = ''.join([merged_quality, reverse_quality[i]])
                continue

            #Error out otherwise
            print("[FASTQ Merge Error] Something went unchecked during the merge.")
            quit()

        return [merged_read, merged_quality]

    def fastq_merge(self):
        """Merge the two files"""

        #Update the user
        print("[FASTQ Merge] Reading in the forward file.")

        #Read the FASTQ files into a list
        print("Opened: " + self.forward_fastq)
        with open(self.forward_fastq, 'r') as forward_file:
            dnaseq = list(islice(forward_file, 1, None, 2))
            
        #Combine the dna sequence and quality together
        list_forward = [[dnaseq[i].rstrip('\n\r'), dnaseq[i+1].rstrip('\n\r')] for i in range(0, len(dnaseq), 2)]

        #Update the user
        print("[FASTQ Merge] Reading in the reverse file.")        

        print("Opened: " + self.reverse_fastq)
        with open(self.reverse_fastq, 'r') as reverse_file:
            dnaseq = list(islice(reverse_file, 1, None, 2))
            
        #Combine the dna sequence and quality together
        list_reverse = [[dnaseq[i].rstrip('\n\r'), dnaseq[i+1].rstrip('\n\r')] for i in range(0, len(dnaseq), 2)]

        #Merge the two lists
        list_fastqs = list(zip(list_forward, list_reverse))

        #Scout 100 reads step 100 to get the mode read start index
        if len(list_fastqs) > 10000:
            max_scouting = 1000#0
        else:
            max_scouting = len(list_fastqs)

        list_scouting = []
        for i in range(0, max_scouting, 100):
            list_scouting.append(list_fastqs[i])

        #Update the user
        print("[FASTQ Merge] Mapping processes to scout 100 test sequences.")

        #Make a pool of worker threads
        pool = ThreadPool(self.processes)

        #Open the DNA sequences in their own threads and return the results
        try:
            scouting_results = pool.map(unwrap_self_merge_pairwise_scout, zip([self]*len(list_scouting), list_scouting))
        except MemoryError:
            print("[FASTQ Merge Error] Out of memory. Please use a 64bit install of python or install more ram.")
            exit(1)

        #Close the pool and make sure the results are in
        pool.close()
        pool.join()

        #Update the user
        print("[FASTQ Merge] Reducing the results to calculate statistics for optimizations.")

        #Create a list of indexes
        list_index = list(x[0] for x in scouting_results if x[1] >= self.min_coverage)

        #Create a list with the mode as the first element then the rest of the values w/in 1 SD
        try:
            index_mode = mode(list_index)
            index_sd = int(stdev(list_index))
            
            #Initialize a new list with the mode as the first entry
            self.list_tight_index = [index_mode]

            #Add the values around the standard devation to the new list
            for i in range(index_mode - index_sd, index_mode + index_sd + 1):
                if i != mode(list_index):
                    self.list_tight_index.append(i)
        except StatisticsError:
            #Fall back to a single starting point
            self.list_tight_index = [1]
        
        #Set the direction mode, 1 = FWD, 2 = REV, 3 = unknown
        try:
            self.direction_mode = mode(list(x[2] for x in scouting_results if x[1] >= self.min_coverage))
        except StatisticsError:
            self.direction_mode = 3

        #Update the user
        print("[FASTQ Merge] Mapping processes to merge reads.")

        #Make a pool of worker threads
        pool = ThreadPool(self.processes)

        #Open the DNA sequences in their own threads and return the results
        try:
            results = pool.map(unwrap_self_merge_pairwise, zip([self]*len(list_fastqs), list_fastqs))
        except MemoryError:
            print("[FASTQ Merge Error] Out of memory. Please use a 64bit install of python or install more ram.")
            exit(1)

        #Close the pool and make sure the results are in
        pool.close()
        pool.join()

        #Update the user
        print("[FASTQ Merge] Reducing the results and writing the merged file.")

        #Write our new file
        file_out = open(self.out_prefix + "_Merge.fastq", 'w')
        counter = 0

        for result in results:
            if result[0] is not None:
                file_out.write(''.join(["@Read: ", str(counter), "\n", result[0], "\n+\n", result[1], "\n"]))
                counter = counter + 1

        file_out.close()

        #Update the user
        print(' '.join(["[FASTQ Merge] Wrote", str(counter), "reads."]))
        
        return "[FASTQ Merge] Wrote " + str(counter) + " reads."

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
    parser = ArgumentParser(description='FASTQ Merge')
    parser.add_argument('-f', dest='forward', action='store', help='Forward file')
    parser.add_argument('-r', dest='reverse', action='store', help='Reverse file')
    parser.add_argument('-j', dest='processes', action='store', help='Number of processes to spawn')
    parser.add_argument('-o', dest='output', action='store', help='Output file path + prefix')
    parser.add_argument('-m', dest='coverage', action='store', help='Minimum read coverage (0.2)')
    args = parser.parse_args()

    merge_settings = {
        'Processes':args.processes,
        'Forward_FASTQ':args.forward,
        'Reverse_FASTQ':args.reverse,
        'Out_Prefix':args.output,
        'Min_Coverage':args.coverage,
        }

    #Create our object and run
    obj_merge = fastq_merge(merge_settings)
    obj_merge.fastq_merge()
