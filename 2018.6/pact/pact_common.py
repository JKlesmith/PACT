#!/usr/bin/env python3
#
# Protein Analysis and Classifier Toolkit
# Author: Justin R. Klesmith
# Copyright (C) 2018 by Regents of the University of Minnesota
# Copyright (C) 2018 by Justin R. Klesmith
#
# This software is released under GNU General Public License 3
# Additional license options available at http://license.umn.edu
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

"""Common functions for the Protein Analysis and Classifier Toolkit"""

from sys import version_info
from math import sqrt, pow

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[PACT Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

def filename_sanitize(filename):
    """Sanitize a string to be used as a filename"""
    keepcharacters = (' ', '.', '_')
    string_file_out = "".join(c for c in filename if c.isalnum() or c in keepcharacters).rstrip()
    return string_file_out

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

def command_line_args(flag):
    """Parse the command line arguments and return the value for a flag"""
    from sys import argv

    #Get the index of our flag then return the value from the next item.
    try:
        return True, argv[argv.index(flag) + 1]
    except ValueError:
        return False, None

def dna_dicts():
    """Return common DNA lookup dicts."""
    DNA_Table = {'TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA',
    'TTG', 'TCG', 'TAG', 'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC',
    'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT',
    'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG',
    'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA',
    'GTG', 'GCG', 'GAG', 'GGG'}

    Translation_Table = {'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
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

    Translation_Table_AtoD = {
    '*' : ['TAA', 'TAG', 'TGA'],
    'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
    'C' : ['TGT', 'TGC'],
    'D' : ['GAT', 'GAC'],
    'E' : ['GAA', 'GAG'],
    'F' : ['TTT', 'TTC'],
    'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
    'H' : ['CAT', 'CAC'],
    'I' : ['ATT', 'ATC', 'ATA'],
    'K' : ['AAA', 'AAG'],
    'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M' : ['ATG'],
    'N' : ['AAT', 'AAC'],
    'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q' : ['CAA', 'CAG'],
    'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
    'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
    'W' : ['TGG'],
    'Y' : ['TAT', 'TAC']}

    AA_Table = '*FWYPMILVAGCSTNQDEHKR'

    return [DNA_Table, Translation_Table, Translation_Table_AtoD, AA_Table]

def pretty_counter_dicts(dict_counter):
    """Convert a counter dict to a pretty line"""
    return ', '.join([class_name + ": " + str(dict_counter[class_name]) 
                         for class_name in sorted(dict_counter)])

def get_bool(input, relation, value):
    """
    Another way to evaluate relations
    Useage: if self.get_bool(int/float_value, str_relation, int/float_value):
    """

    #Setup a relation table
    import operator
    ops = {
    '>': operator.gt,
    '<': operator.lt,
    '>=': operator.ge,
    '<=': operator.le,
    '=': operator.eq,
    '!': operator.ne,
    }
    return ops[relation](input, value)

def create_csv(str_filename, header, text):
    """Create a csv"""

    #Open the output file and write it
    with open(str_filename, "w") as csv_out:

        #Write the header
        csv_out.write(header + "\n")

        #Write the text
        csv_out.write(text + "\n")

    return "Wrote CSV: " + str_filename

def euc_dist(x1, x2, y1, y2, z1, z2):
    """Calculate the Euclidean distance between two points"""
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))

def save_pact_file(data_structure, pact_filename):
    """Save the results of this program into a python friendly form."""

    #Open and write the output file
    from pickle import dump, HIGHEST_PROTOCOL
    with open(pact_filename + '.pact', 'wb') as outfile:
        dump(data_structure, outfile, protocol=HIGHEST_PROTOCOL)

    return "[PACT Output] Output saved as " + pact_filename + '.pact'

def open_pact_file(pact_filename):
    """Load the pact file."""

    from pickle import load
        
    #See if the file exists
    if file_checker(pact_filename + '.pact'):

        #Load the pact file
        with open(pact_filename + '.pact', 'rb') as infile:
            pact_file = load(infile)

    return pact_file