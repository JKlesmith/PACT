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

"""codon_condenser - find the common dengenerate codon for a set of amino acids."""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Codon Condenser Error] Your Python interpreter is too old. Minimium version required is Python 3.4.")
    quit()

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class codon_condenser:
    """Find the common dengenerate codon for a set of amino acids"""

    def __init__(self):
        """Initialize the class varibles"""

        #Dict for the codon to amino acid
        self.translation_table = {
            'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
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

        #Dict for amino acid to codon
        self.dict_amino_to_codon = {
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

        #Dict of the nucleotide codes
        self.dict_ntcodes = {
            'G':['G'], 
            'A':['A'], 
            'T':['T'], 
            'C':['C'],
            'R':['G', 'A'], 
            'Y':['T', 'C'], 
            'M':['A', 'C'], 
            'K':['G', 'T'], 
            'S':['G', 'C'], 
            'W':['A', 'T'], 
            'H':['A', 'C', 'T'], 
            'B':['G', 'T', 'C'], 
            'V':['G', 'C', 'A'], 
            'D':['G', 'A', 'T'], 
            'N':['G', 'A', 'T', 'C']}

        #Dict of DNA base to degenerate base
        self.dict_ntcodes_inv = {
            'G':'G', 'A':'A', 'T':'T', 'C':'C',
            'GA':'R', 'AG':'R',
            'TC':'Y', 'CT':'Y',
            'AC':'M', 'CA':'M',
            'GT':'K', 'TG':'K',
            'GC':'S', 'CG':'S',
            'AT':'W', 'TA':'W',
            'ACT':'H', 'ATC':'H', 'TCA':'H', 'TAC':'H', 'CTA':'H', 'CAT':'H',
            'GTC':'B', 'GCT':'B', 'CGT':'B', 'CTG':'B', 'TCG':'B', 'TGC':'B',
            'GCA':'V', 'GAC':'V', 'CAG':'V', 'CGA':'V', 'AGC':'V', 'ACG':'V',
            'GAT':'D', 'GTA':'D', 'ATG':'D', 'AGT':'D', 'TGA':'D', 'TAG':'D',
            'ACTG':'N', 'ACGT':'N', 'ATGC':'N', 'ATCG':'N', 'AGTC':'N', 'AGCT':'N',
            'CATG':'N', 'CAGT':'N', 'TAGC':'N', 'TACG':'N', 'GATC':'N', 'GACT':'N',
            'CTAG':'N', 'CGAT':'N', 'TGAC':'N', 'TCAG':'N', 'GTAC':'N', 'GCAT':'N',
            'CTGA':'N', 'CGTA':'N', 'TGCA':'N', 'TCGA':'N', 'GTCA':'N', 'GCTA':'N'}

        #String of amino acids codes
        self.aa_table = 'ACDEFGHIKLMNPQRSTVWY*'

        #String of codon codes
        self.codon_codes = 'GATCRYMKSWHBVDN'

        #Dict for the codon usage in yeast
        self.yeast_optimal = {
            "TTT":0.59, "TCT":0.26, "TAT":0.56, "TGT":0.63, 
            "TTC":0.41, "TCC":0.16, "TAC":0.44, "TGC":0.37, 
            "TTA":0.28, "TCA":0.21, "TAA":0.47, "TGA":0.30, 
            "TTG":0.29, "TCG":0.10, "TAG":0.23, "TGG":1.00, 
            "CTT":0.13, "CCT":0.31, "CAT":0.64, "CGT":0.14, 
            "CTC":0.06, "CCC":0.15, "CAC":0.36, "CGC":0.06, 
            "CTA":0.14, "CCA":0.42, "CAA":0.69, "CGA":0.07, 
            "CTG":0.11, "CCG":0.12, "CAG":0.31, "CGG":0.04, 
            "ATT":0.46, "ACT":0.35, "AAT":0.59, "AGT":0.16, 
            "ATC":0.26, "ACC":0.22, "AAC":0.41, "AGC":0.11, 
            "ATA":0.27, "ACA":0.30, "AAA":0.58, "AGA":0.48, 
            "ATG":1.00, "ACG":0.14, "AAG":0.42, "AGG":0.21, 
            "GTT":0.39, "GCT":0.38, "GAT":0.65, "GGT":0.47, 
            "GTC":0.21, "GCC":0.22, "GAC":0.35, "GGC":0.19, 
            "GTA":0.21, "GCA":0.29, "GAA":0.70, "GGA":0.22, 
            "GTG":0.19, "GCG":0.11, "GAG":0.30, "GGG":0.12}
        
        #Dict for the codon usage in human cells
        self.human_optimal = {
            "TTT":0.46,"TCT":0.19,"TAT":0.44,"TGT":0.46,
            "TTC":0.54,"TCC":0.22,"TAC":0.56,"TGC":0.54,
            "TTA":0.08,"TCA":0.15,"TAA":0.30,"TGA":0.47,
            "TTG":0.13,"TCG":0.05,"TAG":0.24,"TGG":1.00,
            "CTT":0.13,"CCT":0.29,"CAT":0.42,"CGT":0.08,
            "CTC":0.20,"CCC":0.32,"CAC":0.58,"CGC":0.18,
            "CTA":0.07,"CCA":0.28,"CAA":0.27,"CGA":0.11,
            "CTG":0.40,"CCG":0.11,"CAG":0.73,"CGG":0.20,
            "ATT":0.36,"ACT":0.25,"AAT":0.47,"AGT":0.15,
            "ATC":0.47,"ACC":0.36,"AAC":0.53,"AGC":0.24,
            "ATA":0.17,"ACA":0.28,"AAA":0.43,"AGA":0.21,
            "ATG":1.00,"ACG":0.11,"AAG":0.57,"AGG":0.21,
            "GTT":0.18,"GCT":0.27,"GAT":0.46,"GGT":0.16,
            "GTC":0.24,"GCC":0.40,"GAC":0.54,"GGC":0.34,
            "GTA":0.12,"GCA":0.23,"GAA":0.42,"GGA":0.25,
            "GTG":0.46,"GCG":0.11,"GAG":0.58,"GGG":0.25}

        return

    def codon_to_wanted(self, codon, amino_acids):
        """Analyse each codon then return the results (called by codon_condense)."""

        #Make the list of possible codons
        #[('G', 'A', 'G'), ('G', 'A', 'C'), ...
        #Then join them 'GAG'...
        list_possibilities = list("".join((x, y, z))
                                  for x in self.dict_ntcodes[codon[0]] 
                                  for y in self.dict_ntcodes[codon[1]] 
                                  for z in self.dict_ntcodes[codon[2]])

        #Let's build a dict of all the possible codons mapped to AA's and stops
        #{'A': [], 'C': [], 'D': [],
        dict_amino = {amino: [] for amino in self.aa_table}

        #Assign the codons to the aminos
        #{'A': ['GCG', 'GCT'], 'C': [], 'D': ['GAT'], 'E': ['GAG'],
        for codon_str in list_possibilities:
            dict_amino[self.translation_table[codon_str]].append(codon_str)

        #Setup the initial counts for variabless
        count_total_aminos = 0
        count_wanted_codons = 0
        str_wanted_codons = ""
        str_all_aminos = ""
        bool_filtered = False

        #Loop the aa dict for total number of aas and make a string of all hits
        for amino in self.aa_table:
            if len(dict_amino[amino]) != 0:
                #Make a string for the total aa hits
                str_all_aminos = str_all_aminos + amino + ":" + str(len(dict_amino[amino])) + ", "
        
                #Count the total number of amino acids
                if amino != "*":
                    count_total_aminos = count_total_aminos + 1

        #Sort alphabetically the input AA's
        sorted_amino_acids = ''.join(sorted(amino_acids))

        #Calculate the % of wanted AA's, and filter our results if we have a codon that doesn't have our wanted
        for i in range(0, len(amino_acids)):
            #Calculate the number of codons per hit
            count_wanted_hits = len(dict_amino[sorted_amino_acids[i]])
        
            #Create a string with the counts of our hits
            str_wanted_codons = str_wanted_codons + sorted_amino_acids[i] + ":" + str(count_wanted_hits) + ", "
    
            #Count the total number of codons that our wanted have
            count_wanted_codons = count_wanted_codons + count_wanted_hits
        
            #If no codons then filter this result out
            if count_wanted_hits == 0:
                bool_filtered = True
    
        #Finish up the math on the %wanted, total number of AA and Stops
        count_total = len(list_possibilities)
        percent_stop = (len(dict_amino['*'])/count_total)*100
        percent_wanted = (count_wanted_codons/count_total)*100

        return (count_total, count_total_aminos, percent_wanted, percent_stop, 
                str_wanted_codons.rstrip(", "), str_all_aminos.rstrip(", "), bool_filtered)
    
    def codon_analysis(self, codon):
        """Input a codon and output the amino acids encoded."""

        #Convert our codon to upper
        codon = codon.upper()

        #Check to see if our string is not three characters
        if len(codon) != 3:
            print("[Codon Condenser Error] The codon input is not three characters.")
            quit()

        #Check to see if the letters are actual codon codes
        for letter in codon:
            if letter not in self.codon_codes:
                print("[Codon Condenser Error] Unknown codon letter: "+letter)
                quit()

        #Now make the list of possible codons
        list_possibilities = list("".join([x, y, z]) for x in self.dict_ntcodes[codon[0]] 
                                  for y in self.dict_ntcodes[codon[1]] 
                                  for z in self.dict_ntcodes[codon[2]])

        #Let's build a dict of amino acid to possible codons
        dict_amino = {amino: [] for amino in self.aa_table}

        for listcodon in list_possibilities:
            dict_amino[self.translation_table[listcodon]].append(listcodon) #'A': ['GCG', 'GCA', 'GCT', 'GCC'], ...

        #Calculate the total number of AA and Stops
        count_codons = len(list_possibilities)
        count_stops = len(dict_amino['*'])
        count_aminos = count_codons - count_stops

        #Print the list of bases at each site from the input codon
        output_str = (
        "[Codon Condenser] Codon to analyze: " + codon + "\n"
        "[Codon Condenser] 1st base: " + "".join(self.dict_ntcodes[codon[0]]) + "\n"
        "[Codon Condenser] 2nd base: " + "".join(self.dict_ntcodes[codon[1]]) + "\n"
        "[Codon Condenser] 3rd base: " + "".join(self.dict_ntcodes[codon[2]]) + "\n"
        "[Codon Condenser] Total number of codons: " + str(count_codons) + "\n"
        "[Codon Condenser] Total number of amino acid encoding codons: " + str(count_aminos) + " ({0:.1f}".format((count_aminos/count_codons)*100) + "%)\n"
        "[Codon Condenser] Total number of stop codons: " + str(count_stops) + " ({0:.1f}".format((count_stops/count_codons)*100) + "%)\n"
        "AA\t#\t%\tCodons\n"
        )
   
        #Let's output the AA to Codon mappings
        for amino in self.aa_table:
            output_str = output_str + "\t".join([amino, 
                                str(len(dict_amino[amino])), 
                                "{0:.1f}".format((len(dict_amino[amino])/count_codons)*100),
                                ",".join(dict_amino[amino])]) + "\n"
        
        #Return the amino acid lookup dict
        return dict_amino, output_str

    def codon_condense(self, amino_acids):
        """Method to find the degenerate codons for input amino acids."""

        #Import the product function
        from itertools import product

        #Convert the input to upper
        amino_acids = amino_acids.upper()

        #Inform the user of their decision    
        output_str = "[Codon Condenser] Amino acids to condense: " + ''.join(sorted(amino_acids)) + "\n"

        #Loop through each amino acid (Ex: for ADS)
        #[['GCT', 'GCC', 'GCA', 'GCG'], ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], ['GAT', 'GAC']]
        try:
            list_codon_groups = list(self.dict_amino_to_codon[amino_acids[i]] for i in range(0, len(amino_acids)))
        except KeyError:
            print("[Codon Condenser Error] Non-amino acid in input.")
            quit()

        #Produce a list of the unique bases for each codon base position seperated by amino acid
        #[['G'], ['T', 'A'], ['G']] for base 0
        codon1 = list(list(set([codon[0] for codon in group])) for group in list_codon_groups)
        codon2 = list(list(set([codon[1] for codon in group])) for group in list_codon_groups)
        codon3 = list(list(set([codon[2] for codon in group])) for group in list_codon_groups)

        #Use the product iterator to make the combinations from the lists of lists
        #[('G', 'A', 'G'), ('G', 'T', 'G')] aka: product(*codon1) * = splat = split the list into the argument field
        #Then join the results as a string and remove duplicates
        #['GT', 'GA']

        codon1_combined = list("".join(set((x))) for x in list(product(*codon1)))
        codon2_combined = list("".join(set((x))) for x in list(product(*codon2)))
        codon3_combined = list("".join(set((x))) for x in list(product(*codon3)))
        
        #Convert the bases to codon possibilities
        #['K', 'R']
        codon1_degcodes = list(set(list(self.dict_ntcodes_inv[base] for base in codon1_combined)))
        codon2_degcodes = list(set(list(self.dict_ntcodes_inv[base] for base in codon2_combined)))
        codon3_degcodes = list(set(list(self.dict_ntcodes_inv[base] for base in codon3_combined)))

        #Test for N
        if len(codon1_degcodes) == 14:
            codon1_degcodes.append('N')

        if len(codon2_degcodes) == 14:
            codon2_degcodes.append('N')
        
        if len(codon3_degcodes) == 14:
            codon3_degcodes.append('N')
        
        output_str = output_str + "[Codon Condenser] 1st base option(s): " + "".join(codon1_degcodes) + "\n"
        output_str = output_str + "[Codon Condenser] 2nd base option(s): " + "".join(codon2_degcodes) + "\n"
        output_str = output_str + "[Codon Condenser] 3rd base option(s): " + "".join(codon3_degcodes) + "\n"

        #From the codes create different combinations across all three positions
        #[('R', 'M', 'H'), ('R', 'M', 'S'), ('R', 'M', 'T'), ('R', 'M', 'W'), ('R', 'M', 'Y') ...
        #Then join them ['RVT', 'RVW', 'RVY',...
        list_whole_codon = list("".join((x, y, z)) for x in codon1_degcodes for y in codon2_degcodes for z in codon3_degcodes)

        #Assess each codon
        #{'RMM': ['RMM', 8, 6, 37.5, 0.0, 'A:2, D:1, S:0', 'A:2, D:1, E:1, K:1, N:1, T:2', True], ...
        #[UniqueCodon] = Codon, TotalCodons, TotalAAs, %AAWanted, %Stop, #CodonsWanted, All AA and # of each, Filter Out this result
        dict_codon_analysis = {combo: [combo] + list(self.codon_to_wanted(combo, amino_acids)) for combo in list_whole_codon}

        #Sort our results
        #[{'codon': 'RVY', 'totalc': 12, 'totalaa': 6, 'pcwanted': 50.0, 'pstop': 0.0, 'ncwanted':...
        list_output = list({'codon':dict_codon_analysis[dictcodon][0], 
                            'totalc':dict_codon_analysis[dictcodon][1], 
                            'totalaa':dict_codon_analysis[dictcodon][2], 
                            'pcwanted':dict_codon_analysis[dictcodon][3], 
                            'pstop':dict_codon_analysis[dictcodon][4], 
                            'ncwanted':dict_codon_analysis[dictcodon][5], 
                            'allaa':dict_codon_analysis[dictcodon][6]}
                            for dictcodon in dict_codon_analysis 
                            if dict_codon_analysis[dictcodon][7] is False)

        #Sort our table on the %wanted, totalAA and %stop
        list_output.sort(key=lambda k : (-k['pcwanted'], k['totalaa'], k['pstop']))

        #Re-map into a new list object with the header in the first row
        list_pretty_output = []
        if len(amino_acids) is 1:
            list_pretty_output.append(["Codon", 
                                "Total Codons", 
                                "Total # AA", 
                                "% Codons of Wanted", 
                                "% Stop Codons", 
                                "# Codons of Wanted", 
                                "# Codons of All", 
                                "Yeast Usage", 
                                "Human Usage"])
            for outitem in list_output:
                list_pretty_output.append([outitem['codon'], 
                                    str(outitem['totalc']), 
                                    str(outitem['totalaa']), 
                                    "{0:.1f}".format(outitem['pcwanted']), 
                                    "{0:.1f}".format(outitem['pstop']), 
                                    outitem['ncwanted'], 
                                    outitem['allaa'], 
                                    "{0:.2f}".format(self.yeast_optimal[outitem['codon']]), 
                                    "{0:.2f}".format(self.human_optimal[outitem['codon']])])      
        else:
            list_pretty_output.append(["Codon", 
                                "Total Codons", 
                                "Total # AA", 
                                "% Codons of Wanted", 
                                "% Stop Codons", 
                                "# Codons of Wanted", 
                                "# Codons of All"])
            for outitem in list_output:
                list_pretty_output.append([outitem['codon'], 
                                    str(outitem['totalc']), 
                                    str(outitem['totalaa']), 
                                    "{0:.1f}".format(outitem['pcwanted']), 
                                    "{0:.1f}".format(outitem['pstop']), 
                                    outitem['ncwanted'], 
                                    outitem['allaa']])  
     
        #Figure out the max length of each column (inital the vars with len of 0)
        column_len = {j: 0 for j in range(0, len(list_pretty_output[0]))}

        #Figure out the max length of each column
        for row in list_pretty_output: #Iter each row
            for i in range(0, len(row)): #Iter each column
                if len(row[i]) > column_len[i]:
                    column_len[i] = len(row[i])

        #Print out our results (adjusting for the column width)
        for row in list_pretty_output:
            output = ""
            for i in range(0, len(row)): #Iter each column
                output = output + "".join(row[i].ljust(column_len[i])) + "\t"
            output_str = output_str + output + "\n"
        
        return dict_codon_analysis, output_str

if __name__ == '__main__':
    from argparse import ArgumentParser
    
    #Parse the command line for inputs
    parser = ArgumentParser(description='CodonCondenser')
    parser.add_argument('-a', dest='aminoacid', action='store', help='String of AAs to encode')
    parser.add_argument('-c', dest='codon', action='store', help='Codon to calculate AA usage')
    args = parser.parse_args()

    #Print a header
    print("[Codon Condenser] Codon Condenser\n"
            "[Codon Condenser] Author: "+__author__+"\n"
            "[Codon Condenser] Contact: "+__email__[0]+", "+__email__[1]+"\n"
            "[Codon Condenser] "+__copyright__+"\n"
            "[Codon Condenser] Version: "+__version__+"\n"
            "[Codon Condenser] License: "+__license__+"\n"
            "[Codon Condenser] Github [user: JKlesmith] (www.github.com)")

    #Create our object
    ccobj = codon_condenser()
    
    #Check to see if we have a amino acid or a codon
    if args.aminoacid is None and args.codon is None:
        print("[Codon Condenser Error] Missing flags.")
        quit()
          
    #If we have a amino acid defined and make it uppercase
    if args.aminoacid is not None and args.aminoacid != "":
        dict_aminos, output_str = ccobj.codon_condense(args.aminoacid)
        print(output_str)

    #If we have a codon defined and make it uppercase
    if args.codon is not None and args.codon != "":
        dict_codon, output_str = ccobj.codon_analysis(args.codon)
        print(output_str)
