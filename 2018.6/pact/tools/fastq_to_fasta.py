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

"""fastq_to_fasta - convert a fastq file to a fasta file format."""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Fastq to Fasta Error] Your Python interpreter is too old. Minimium version required is Python 3.4.")
    quit()

from itertools import islice

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class fastq_to_fasta():
    """Convert a fastq file to a fasta file"""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""
        
        #Check to see if the file exists
        if self.file_checker(settings_dict["fastq_file"]):
            self.file_fastq = settings_dict["fastq_file"]

        return

    def file_checker(file_location):
        """Check to see if the file exists"""
        from os.path import isfile

        try:
            if isfile(file_location):
                return True
            else:
                print("[File Import Error] The file: "+file_location+" cannot be found.")
                quit()
                return False
        except TypeError:
                print("[File Import Error] The file: "+str(file_location)+" cannot be found.")
                quit()
                return False

    def fastq_to_fasta(file_fastq):
        """Convert the fastq file to a fasta file per read"""

        dna_list = []
        #Open the input file into a list
        with open(self.file_fastq, 'r') as file_in:
            dnaseq = islice(file_in, 1, None, 4)
            for line in dnaseq:
                dna_list.append(line)

        #Write the AA mappings to a file
        counter = 0
        with open(self.file_fastq + ".fasta", 'w') as output_fastq:
            for sequence in dna_list:
                output_fastq.write(">" + str(counter) + "\n")
                output_fastq.write(sequence + "\n")
                counter = counter + 1

        return

if __name__ == '__main__':
    from argparse import ArgumentParser

    #Parse the command line for inputs
    parser = ArgumentParser(description='FASTQ to FASTA')
    parser.add_argument('-f', dest='filename', action='store', required=True, help='Run Mode')
    args = parser.parse_args()

    #Print a header
    print("[Fastq to Fasta] Fastq to Fasta\n"
            "[Fastq to Fasta] Author: "+__author__+"\n"
            "[Fastq to Fasta] Contact: "+__email__[0]+", "+__email__[1]+"\n"
            "[Fastq to Fasta] "+__copyright__+"\n"
            "[Fastq to Fasta] Version: "+__version__+"\n"
            "[Fastq to Fasta] License: "+__license__+"\n"
            "[Fastq to Fasta] Github [user: JKlesmith] (www.github.com)")

    obj_ftf = fastq_to_fasta({"fastq_file":args.filename})
    obj_ftf.fastq_to_fasta()
