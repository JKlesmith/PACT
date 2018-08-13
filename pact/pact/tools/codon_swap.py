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

"""codon_swap - swap the codons to synonymous codons"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Primer Design Error] Your Python interpreter is too old. Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class codon_swap:
    """The synonymous codon swap class."""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""

        self.aa_table = '*ACDEFGHIKLMNPQRSTVWY'
        self.codon_table = {'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
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
        self.ecoli_table = {'TTT':'TTC', 'TCT':'AGC', 'TAT':'TAC', 'TGT':'TGC',
        'TTC':'TTT', 'TCC':'AGC', 'TAC':'TAT', 'TGC':'TGT',
        'TTA':'CTG', 'TCA':'AGC', 'TAA':'TAA', 'TGA':'TAA',
        'TTG':'CTG', 'TCG':'AGC', 'TAG':'TGA', 'TGG':'TGG',
        'CTT':'CTG', 'CCT':'CCG', 'CAT':'CAC', 'CGT':'CGC',
        'CTC':'CTG', 'CCC':'CCG', 'CAC':'CAT', 'CGC':'CGT',
        'CTA':'CTG', 'CCA':'CCG', 'CAA':'CAG', 'CGA':'CGT',
        'CTG':'CTG', 'CCG':'CCA', 'CAG':'CAA', 'CGG':'CGC',
        'ATT':'ATC', 'ACT':'ACC', 'AAT':'AAC', 'AGT':'AGC',
        'ATC':'ATT', 'ACC':'ACG', 'AAC':'AAT', 'AGC':'TCT',
        'ATA':'ATT', 'ACA':'ACC', 'AAA':'AAA', 'AGA':'CGT',
        'ATG':'ATG', 'ACG':'ACC', 'AAG':'AAA', 'AGG':'CGC',
        'GTT':'GTG', 'GCT':'GCG', 'GAT':'GAC', 'GGT':'GGC',
        'GTC':'GTA', 'GCC':'GCG', 'GAC':'GAT', 'GGC':'GGT',
        'GTA':'GTC', 'GCA':'GCG', 'GAA':'GAG', 'GGA':'GGC',
        'GTG':'GTT', 'GCG':'GCC', 'GAG':'GAA', 'GGG':'GGT'}

        self.wildtype = settings_dict['WT_DNA'].rstrip('\r\n').upper()
        #self.organism = settings_dict['Organism'].rstrip('\r\n').lower()
        return

    def codon_swap(self):
        wildtype_aa_seq = ""
        new_dna_seq = ""
        new_aa_seq = ""

        len_wildtype = int(len(self.wildtype)/3)
        if len_wildtype % 3 != 0:
            print("DNA sequence is not a mulitple of three.")
            exit(1)

        for i in range(0, len_wildtype):
            wildtype_aa_seq = wildtype_aa_seq + self.codon_table[self.wildtype[i*3:i*3+3]]
            new_dna_seq = new_dna_seq + self.ecoli_table[self.wildtype[i*3:i*3+3]]
            new_aa_seq = new_aa_seq + self.codon_table[self.ecoli_table[self.wildtype[i*3:i*3+3]]]

        string_output = ("Wild-type DNA sequence\n"
        ""+self.wildtype+"\n"
        "Codon swapped DNA sequence\n"
        ""+new_dna_seq+"\n"
        "Wild-type amino acid sequence\n"
        ""+wildtype_aa_seq+"\n"
        "Codon swapped amino acid sequence\n"
        ""+new_aa_seq)

        return new_dna_seq, string_output

if __name__ == '__main__':
    
    #Get commandline arguments
    from argparse import ArgumentParser
    
    parser = ArgumentParser(description='CodonSwap: for swapping codons to synomymous mutations')
    #parser.add_argument('-o', dest='organism', action='store', nargs='?', const=1, default='', help='Organism type.')
    parser.add_argument('-w', dest='wildtype', action='store', nargs='?', const=1, default='', help='Wild-type DNA sequence.')
    args = parser.parse_args()

    #Print a header
    print("[Codon Swap] Codon Swap\n"
            "[Codon Swap] Author: "+__author__+"\n"
            "[Codon Swap] Contact: "+__email__[0]+", "+__email__[1]+"\n"
            "[Codon Swap] "+__copyright__+"\n"
            "[Codon Swap] Version: "+__version__+"\n"
            "[Codon Swap] License: "+__license__+"\n"
            "[Codon Swap] Github [user: JKlesmith] (www.github.com)")

    #Create our object
    cs = codon_swap({"WT_DNA":args.wildtype})
    new_dna_seq, string_output = cs.codon_swap()

    print(string_output)
