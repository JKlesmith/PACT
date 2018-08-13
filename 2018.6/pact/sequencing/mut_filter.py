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

"""Mutation Filter - filter out mutations"""

from sys import version_info, maxsize

#Check to see if we're on 64-bit
if maxsize < 2**32:
    print("[PACT WARNING] This is not a 64-bit system, memory issues may occur with large files!")

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Mutation Filter Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from multiprocessing import Pool as ThreadPool
from multiprocessing import cpu_count
from pact.pact_common import file_checker

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

#Multiprocessing doesn't work inside of classes (workaround to expose the method as a top level function)
def unwrap_self_design_line_filter(arg, **kwarg):
    return mut_filter.design_line_filter(*arg, **kwarg)

class mut_filter:
    """Fiter out reads that have mutations outside of our library."""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""

        #Assign our varibles to the class
        self.wtaa = settings_dict['WTAA']
        self.FirstAAMutated = int(settings_dict['FirstAAMutated'])-1

        #Check to see if the number of processes is logical
        if settings_dict['Processes'] <= 0 or settings_dict['Processes'] is None:
            self.processes = cpu_count()
        else:
            self.processes = settings_dict['Processes']

        #Check to see if the input file exists
        if file_checker(settings_dict['In_File']):
            self.in_file = settings_dict['In_File']

        #Setup the output prefix
        self.out_prefix = settings_dict['Out_Prefix']

        #Check to see if the mut count threshold is logical
        if int(settings_dict['MutThreshold']) < 1:
            self.mutcountthreshold = 1
        else:
            self.mutcountthreshold = int(settings_dict['MutThreshold'])

        #Import the mutation design list
        self.list_mutation_design = []
        for group in literal_eval(settings_dict['mutcodons']):
            #If we have n between two numbers
            if 'n' in group:
                #We need exactly three to expand the range
                if len(group) == 3:
                    #Get the start and end position, then iterate through the two points
                    self.list_mutation_design.append(list(i for i in range(group[0], group[-1] + 1)))
                else:
                    print("[Error] The codon list is missing either a start or end point around the 'n' location.")
                    exit(1)
            else:
                #Add the other ranges without changing them
                self.list_mutation_design.append(group)
        return

    def design_line_filter(self, line):
        """Filter each sequence line to see if it is in our designs"""
        
        #Split the tab
        splitline = line.split('\t')

        #Return wild-type sequences
        if splitline[0] == "WT":
            return ['\t'.join(["[WT]", splitline[3]]), True, False, False]

        list_locations = list(map(int, splitline[0].split(','))) #Get the locations
        list_mutations = splitline[1].split(',') #Get the mutations
        count_numberofmuts = len(list_locations)

        #Populate the group matches list
        list_group_matches = [None] * count_numberofmuts

        #Match the mutation to a group in mutation location list
        for j in range(0, count_numberofmuts):
            for group in self.list_mutation_design:
                for location in group:
                    if location == list_locations[j]+self.FirstAAMutated+1:
                        list_group_matches[j] = group

        #Check to see if all of the mutations groups match each other
        #Single mutations will pass the above check so check for none
        if (list_group_matches.count(list_group_matches[0]) != count_numberofmuts or
                list_group_matches[0] is None or int(splitline[4].rstrip('\n')) > self.mutcountthreshold):
            return ['\t'.join([str(list_group_matches), splitline[0], splitline[1]]), False, True, False]

        #Create a lookup table based on the designed codons and populate with WT
        matrix_mutations = {}
        string_wt_aa = ""
        for location in list_group_matches[0]:
            matrix_mutations[location-1] = self.wtaa[location-1]
            string_wt_aa = ','.join([string_wt_aa, self.wtaa[location-1]])

        #Parse the string and assign the locations and mutations into the table
        for i in range(0, count_numberofmuts):
            matrix_mutations[list_locations[i]+self.FirstAAMutated] = list_mutations[i]

        #Build the optput strings
        locationoutput = ""
        mutsoutput = ""
        for location in list_group_matches[0]:
            locationoutput = ','.join([locationoutput, str(location-1)])
            mutsoutput = ','.join([mutsoutput, matrix_mutations[location-1]])

        output_filtered = ('\t'.join([str(list_group_matches[0]), locationoutput.lstrip(","), string_wt_aa.lstrip(","),
            mutsoutput.lstrip(","), splitline[0], splitline[1].rstrip(","), splitline[2], splitline[3]]))

        return [output_filtered, False, False, True]

    def mut_filter(self):
        """Filter out the reads that are not in our designed regions."""

        #Update the user
        print("[Mutation Filter] Reading in the FASTQ file.")

        #Open the read FASTQ file
        with open(self.in_file, 'r') as infile:
            list_lines = infile.read().splitlines()[1:]

        #Update the user
        print("[Mutation Filter] Mapping processes to filter the reads.")

        #Make a pool of worker threads
        pool = ThreadPool(self.processes)

        #Open the DNA sequences in their own threads and return the results
        try:
            results = pool.map(unwrap_self_design_line_filter, zip([self]*len(list_lines), list_lines))
        except MemoryError:
            print("[Error] Out of memory. Please use a 64bit install of python or install more ram.")

        #Close the pool and make sure the results are in and mark the list_lises to be reused       
        pool.close()
        pool.join()

        #Update the user
        print("[Mutation Filter] Reducing the results and writing the filtered TSV files.")

        #Open the files
        file_filter_wildtype = open(self.out_prefix + "_Filtered_WildType.tsv", 'w')
        file_filter_rejected = open(self.out_prefix + "_Filtered_Rejected.tsv", 'w')
        file_filter = open(self.out_prefix+"_Filtered.tsv", 'w')

        #Setup the variables for the output
        count_rejected = 0 #Count of mutations outside of our regions
        count_accepted = 0 #Count of correct groups
        count_wildtype = 0 #Count of wild-type sequences

        file_filter_rejected.write("ID\tLocation\tMutations\n")
        file_filter.write("ID\tFull_Location\tWild_Type\tFull_Mutations\tLocation\tMutations\tAA_Sequence\tDNA_Sequence\n")
        file_filter_wildtype.write("ID\tDNA_Sequence\n")

        #Loop each result
        for result in results:
            if result[1]:
                #Wild-type
                count_wildtype = count_wildtype + 1
                file_filter_wildtype.write(''.join([result[0], "\n"]))

            elif result[2]:
                #Rejected Read
                count_rejected = count_rejected + 1
                file_filter_rejected.write(''.join([result[0], "\n"]))

            elif result[3]:
                #Accepted Read
                count_accepted = count_accepted + 1
                file_filter.write(''.join([result[0], "\n"]))

        #Close the file
        file_filter_rejected.close()
        file_filter.close()
        file_filter_wildtype.close()

        output_string = ("[Mutation Filter] Output file prefix: " + self.out_prefix + "\n"
                        "[Mutation Filter] Sequences rejected: " + str(count_rejected) + "\n"
                        "[Mutation Filter] Sequences accepted: " + str(count_accepted) + "\n"
                        "[Mutation Filter] Wild-Type sequences: " + str(count_wildtype) + "\n")
        print(output_string)

        return output_string

if __name__ == '__main__':
    
    #Remind the user that the module needs to be ran within the context of PACT
    print("[Mutation Filter Error] This module needs to be ran within the context of PACT.")
