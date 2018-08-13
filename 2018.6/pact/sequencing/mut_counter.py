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

"""Mutation counter - count filtered mutations"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Mutation Counter Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from collections import Counter
from errno import EINVAL
from pact.pact_common import file_checker

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class mut_counter:
    """Count our filtered library reads."""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""

        #Setup the output prefix
        self.out_prefix = settings_dict['Out_Prefix']

        return

    def mutation_counter(self, file_in, file_out, header):
        """Open our filtered file and then count our mutations"""

        #Check if the file exists
        if file_checker(file_in) == False:
            print("[Mutation Counter Error] Filtered tsv is missing and cannot be counted.")
            quit()

        #Open the input file and count each unique line
        try:
            #Setup a string to keep the remainder in
            int_fileend = 1
            int_header = 1
            str_remainder = ''
            list_lines = []

            #Open our file and read in chunks
            with open(file_in, 'r') as infile:

                #Old Method
                #list_lines = infile.read().splitlines()[1:]
                
                while int_fileend:
                    #Prepend our remainder and read the next chunk
                    block = str_remainder + infile.read(100000)

                    #EOF
                    if len(block) < 100000:
                        int_fileend = 0

                    #Split our line
                    splitline = block.split("\n")

                    #Add our line to the main list (minus the last line if we're not done)
                    for i in range(int_header, len(splitline) - int_fileend):

                        #Set the header to 0 (as we want to skip the first line with it 1 originally)
                        int_header = 0

                        #Append and remove
                        list_lines.append(splitline[i])

                    #Add our remainder
                    str_remainder = '\n'.join(splitline[-1:])

            #Count our lines, Remove empty lines
            counts = Counter(l for l in list(filter(None, list_lines)))

            #Output the existing line and append the count to the end
            with open(file_out, 'w') as file_counted:
                file_counted.write(header)
            
                for line, count in counts.most_common():
                    file_counted.write(line + "\t" + str(count) + "\n")

        except MemoryError:
            print("[Mutation Counter Error] Out of memory. Please install more ram or use 64-bit python.")
            quit()

        except OSError as error:
            #OSError 22 - Invalid Argument, when opening huge files esp. on OSX, will crash on rejected sequences
            if error.errno == EINVAL:
                print("[Mutation Counter Error] Are you on OSX? Bug with Python in opening large files.")
                quit()

        return

    def mut_counter(self):
        """Count our mutations"""

        #Update the user
        print("[Mutation Counter] Counting the mutations.")

        #Count the file with the accepted mutations
        self.mutation_counter(self.out_prefix + "_Filtered.tsv", self.out_prefix + "_Counted.tsv", 
                              "ID\tFull_Location\tWild_Type\tFull_Mutations\tLocation\tMutations\tAA_Sequence\tDNA_Sequence\tCounts\n")

        #Count the file with the rejected mutations
        self.mutation_counter(self.out_prefix + "_Filtered_Rejected.tsv", self.out_prefix + "_Counted_Rejected.tsv", 
                              "ID\tLocation\tMutations\tCounts\n")

        #Count the wild-type synonymous mutation file
        self.mutation_counter(self.out_prefix + "_Filtered_WildType.tsv", self.out_prefix + "_Counted_WildType.tsv",
                              "ID\tDNA_Sequence\tCounts\n")

        #Construct the output string
        output_string = ("[Mutation Counter] Counted the file with the accepted mutations\n"
                        "[Mutation Counter] Counted the file with the rejected mutations\n"
                        "[Mutation Counter] Counted the file with the synonymous wild-type mutations")

        #Update the user
        print(output_string)

        return output_string

if __name__ == '__main__':
    
    #Remind the user that the module needs to be ran within the context of PACT
    print("[Mutation Counter Error] This module needs to be ran within the context of PACT.")
