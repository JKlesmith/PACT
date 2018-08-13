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

"""fastq_split.py - Split a fastq file based on the first"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Fastq Split Error] Your Python interpreter is too old. Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

#Other Imports
from os import stat
from contextlib import ExitStack

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class fastq_split:

    def __init__(self, settings_dict):
        """Initialize the class varibles"""

        #Set the input files
        if self.file_checker(settings_dict['forward_file']):
            self.forward_file = settings_dict['forward_file']

        if self.file_checker(settings_dict['reverse_file']):
            self.reverse_file = settings_dict['reverse_file']

        #Check to see if the two files are the same size        
        if stat(self.forward_file).st_size != stat(self.reverse_file).st_size:
            print('[Fastq Split Error] File sizes unequal.')
            quit()

        #Set the output directory
        self.output_prefix = settings_dict['output_prefix']

        #Cutoff is 400 reads to save blowing up the directory
        self.cutoff = 400

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

    def fastq_split(self):
        """Split the two files"""

        #Open the two files and load the reads into lists
        with open(self.forward_file, 'r') as infile_fwd:
            lines_fwd = infile_fwd.read().splitlines()

        with open(self.reverse_file, 'r') as infile_rev:
            lines_rev = infile_rev.read().splitlines()

        #Calculate the total number of reads
        count_reads = int(len(lines_fwd) / 4)

        #Preprocess the list to find 10-mers that we should store
        dict_10mer = {}
        for i in range(0, count_reads):

            #Get the first 10 bases of the forward read
            mer = lines_fwd[1 + (i * 4)][:10]

            #Add to the mer-key count
            if mer in dict_10mer:
                dict_10mer[mer]["COUNT"] = dict_10mer[mer]["COUNT"] + 1
            else:
                dict_10mer[mer] = {}
                dict_10mer[mer]["COUNT"] = 1

        #Fiter the list for only the 10mer that we should store
        for key in list(dict_10mer.keys()):

            #Delete the key if below the cutoff
            if dict_10mer[key]["COUNT"] < self.cutoff:
                del dict_10mer[key]

        #Open the file handles
        with ExitStack() as cm:
            #Open the forward and reverse files
            for key in dict_10mer:
                dict_10mer[key]["FWD"] = cm.enter_context(open(self.output_prefix + key + '-FWD.fastq', 'a'))
                dict_10mer[key]["REV"] = cm.enter_context(open(self.output_prefix + key + '-REV.fastq', 'a'))

            #Loop through and append to the correct file
            for i in range(0, count_reads):             
                #Get the first 10 bases of the forward read
                key = lines_fwd[1 + (i * 4)][:10]

                #Write to our files
                if key in dict_10mer:
                    dict_10mer[key]["FWD"].write(lines_fwd[0 + (i * 4)] + '\n' + lines_fwd[1 + (i * 4)] 
                                                 + '\n+\n' + lines_fwd[3 + (i * 4)] + '\n')
                    dict_10mer[key]["REV"].write(lines_rev[0 + (i * 4)] + '\n' + lines_rev[1 + (i * 4)] 
                                                 + '\n+\n' + lines_rev[3 + (i * 4)] + '\n')

        #Output Notice
        output_str = "[fastq_split] Split: " + str(count_reads) + " reads."
        print(output_str)

        return output_str

if __name__ == '__main__':
    """Entrypoint if ran directly"""

    #Import our argparse
    from argparse import ArgumentParser

    #Setup our commnand line parser
    parser = ArgumentParser(description='fastq_split: split a fastq file into a common start sequence')
    parser.add_argument('-f', dest='forward', action='store', nargs='?', const=1, required=True, help='Forward File')
    parser.add_argument('-r', dest='reverse', action='store', nargs='?', const=1, required=True, help='Reverse File')
    args = parser.parse_args()

    #Print a header
    print("[Fastq Split] Fastq Split\n"
            "[Fastq Split] Author: "+__author__+"\n"
            "[Fastq Split] Contact: "+__email__[0]+", "+__email__[1]+"\n"
            "[Fastq Split] "+__copyright__+"\n"
            "[Fastq Split] Version: "+__version__+"\n"
            "[Fastq Split] License: "+__license__+"\n"
            "[Fastq Split] Github [user: JKlesmith] (www.github.com)")

    #Send our settings
    settings_dict = {
        "forward_file":args.forward,
        "reverse_file":args.reverse,
        "output_prefix":"./",
        }

    #Setup our main class and call the main routine
    obj_split = fastq_split(settings_dict)
    obj_split.fastq_split()
