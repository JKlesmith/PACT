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

"""homology_classifier - sequence homology for protein analysis"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Homology Classifier Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from sys import platform
from pact.pact_common import file_checker, pretty_counter_dicts

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class homology_classifier:

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("blastp_align_filter"):
            print("[Homology PSSM Error] the config file is missing the section [blastp_align_filter]")

        if not self.config_file.has_section("pssm"):
            print("[Homology PSSM Error] the config file is missing the section [pssm]")

        #Classifier specific variables
        self.aa_table = 'FWYPMILVAGCSTNQDEHKR'
        self.wtaa = self.config_file.get("global", "wtaa").upper()
        self.wtlen = len(self.wtaa)
        self.proteinname = self.config_file.get("global", "output_prefix")
        self.directory = self.dict_protocolconfig['directory']

        return

    def xml_to_fasta(self):
        """Convert the XML file into a fasta file for CDHIT"""

        #Assign these but catch value errors
        try:
            minquerylen = float(self.config_file.get("blastp_align_filter", "minquerylen"))
            minseqid = float(self.config_file.get("blastp_align_filter", "minseqid"))
        except ValueError:
            print("Incorrect value for the blastp_align_filter settings (string entered when it should be a number).")
            quit()

        #Get the input files
        self.file_ncbixml = self.config_file.get("blastp_align_filter", "ncbi_xml")

        #Import our xml tools
        from xml.etree import cElementTree

        #Check inputs for valid values
        file_checker(self.directory + self.file_ncbixml)

        if minquerylen > 1:
            print("[xml_to_fasta] The minimum lenth to query length is >1, setting to 1. (Valid values: 0.0 to 1.0)")
            minquerylen = 1
        elif minquerylen < 0:
            print("[xml_to_fasta] The minimum length to query length is <0, setting to 0. (Valid values: 0.0 to 1.0)")
            minquerylen = 0

        if minseqid > 1:
            print("[xml_to_fasta] The minimum sequence identity is >1, setting to 1. (Valid values: 0.0 to 1.0)")
            minseqid = 1
        elif minseqid < 0:
            print("[xml_to_fasta] The minimum sequence identity is <0, setting to 0. (Valid values: 0.0 to 1.0)")
            minseqid = 0
    
        #Flag to see if our WT is in the XML input
        bool_wt_in_ncbi = False 

        #Import NCBI information from TSV
        with open(self.directory + self.proteinname + ".fa", 'w') as file_fa:
            
            #Parse the XML file
            for event, elem in cElementTree.iterparse(self.directory + self.file_ncbixml):
                
                #Parse each hit
                if elem.tag == "Hit":
                    hit_id = elem.find("Hit_id").text
                    hit_name = elem.find("Hit_def").text
                    hit_accession = elem.find("Hit_accession").text
                    hit_identities = int(elem.find("Hit_hsps").find("Hsp").find("Hsp_identity").text)
                    hit_alignlen = int(elem.find("Hit_hsps").find("Hsp").find("Hsp_align-len").text)
                    hit_sequence = elem.find("Hit_hsps").find("Hsp").find("Hsp_hseq").text
            
                    #Verify that the alignment length is >= 60% of query
                    if float(hit_alignlen/self.wtlen) < minquerylen:
                        elem.clear()
                        continue

                    #Verify that the sequence identity is >= 30%
                    if float(hit_identities/hit_alignlen) < minseqid:
                        elem.clear()
                        continue

                    #If there is non-standard amino acids, remove
                    if any(c in hit_sequence for c in ['B', 'J', 'O', 'U', 'X', 'Z']):
                        elem.clear()
                        continue

                    #Remove dashes from the ncbi match (cd-hit errors them out as bad formatting)
                    hit_sequence_nodash = "".join(char for char in hit_sequence if char != "-")

                    #Check to see if the sequence was WT
                    if hit_sequence_nodash == self.wtaa:
                        bool_wt_in_ncbi = True

                        #Write the wild-type to the file
                        file_fa.write(self.proteinname + '\n')
                        file_fa.write(hit_sequence_nodash + '\n\n')
                    else:
                        #Write the other matches to the file
                        file_fa.write(">" + hit_accession + "_" + hit_id + "_" + hit_name + '\n')
                        file_fa.write(hit_sequence_nodash + '\n\n')
            
                    elem.clear()

            #If WT was not in the NCBI input then add it
            if bool_wt_in_ncbi is False:
                file_fa.write(self.proteinname + '\n')
                file_fa.write(self.wtaa + '\n')

        return

    def cdhit_wt_check(self):
        """Checks to see if WT is still there after CDHIT"""

        #Check inputs for valid values
        file_checker(self.directory + self.proteinname + ".afa")

        #Flag to see if our WT is in the input
        bool_wt_in_ncbi = False

        #Open the file and check it line by line
        with open(self.directory + self.proteinname + ".afa", 'r') as infile:
            for line in infile:

                #Check to see if the sequence was WT
                if line.rstrip('\n\r') == self.proteinname:
                    bool_wt_in_ncbi = True

        #If WT was not in the NCBI input then add it
        if bool_wt_in_ncbi == False:
            with open(self.directory + self.proteinname + ".afa", 'a') as outfile:
                outfile.write('>' + self.proteinname + '\r')
                outfile.write('\n')
                outfile.write(self.wtaa + '\r')

        return

    def process_msa(self):
        """Convert the MSA for use with PSI-Blast"""

        #ProcessMSA Settings
        try:
            nummaxhits = int(self.config_file.get("blastp_align_filter", "nummaxhits"))
        except ValueError:
            print("[Homology Error] An incorrect value in the input for max num hits, will be set as unlimited.")
            nummaxhits = -1

        #Check inputs for valid values
        file_checker(self.directory + self.proteinname + ".msa")

        #Step one: Import msa alignment from MUSCLE and make it one line per sequence
        with open(self.directory + self.proteinname + ".msa", 'r') as infile:

            alignment = ""
            output = ""
            for line in infile:
                #Check to see if we have a header
                if line[0] == ">":
                    #Ignores empty output
                    if len(output) > 0:
                        #Add the current output to the growing alignment varible
                        alignment = alignment + output + "\n"

                    #Empty the current alignment
                    output = ""

                    #Assemble the first line of the new sequence
                    output = output + line.rstrip('\n') + "\t"
                else:
                    #Keep assembling the line
                    output = output + line.rstrip('\n')

        #Step two: Import MSA into a lookup table
        dict_msa_table = {line.split("\t")[0]:line.split("\t")[1].rstrip("\n")
                         for line in alignment.split("\n")
                         if len(line) > 10}

        #Step three: Mark the insertions with the letter Z (fixed from X as X is also used in the blast hits)
        #Step four: Delete the insertions
        string_wt_msa = dict_msa_table[">" + self.proteinname] #Get the wild-type MSA sequence (with insertions)

        #Loop the entries in the dict
        for msa in dict_msa_table:

            #Go through the length of the entire msa string and mark it Z
            temp_zstr = "".join("Z" if string_wt_msa[i] == "-" 
                                else dict_msa_table[msa][i] 
                                for i in range(0, len(string_wt_msa)))

            #Now let's delete the Z'ed characters
            temp_zstr = "".join(char for char in temp_zstr if char != "Z")
        
            #Update the dict
            dict_msa_table[msa] = temp_zstr

        #Step five: Put wild-type on top, re-order the sequences by completeness, and cap the hits for psi-blast
        list_msa_ordered = []
        for msa in dict_msa_table:
            
            #Set the WT to the top and then add the rest into a list with their counts of dashes
            if dict_msa_table[msa] == self.wtaa:
                #Move the wild-type to the top
                list_msa_ordered.insert(0, {'sequence' : dict_msa_table[msa], 
                                            'counts' : dict_msa_table[msa].count('-'),
                                            'wt' : True})
            else:
                list_msa_ordered.append({'sequence' : dict_msa_table[msa], 
                                         'counts' : dict_msa_table[msa].count('-'),
                                         'wt' : False})

        #Sort the table into a new list by counts and wt
        list_msa_ordered.sort(key=lambda k : (k['counts'], -k['wt']))

        #Create a list of just the sequences up to our limit
        if nummaxhits == -1 or nummaxhits < 1: #Unlimited number of output sequences
            list_msa_pb = list(x['sequence'] for x in list_msa_ordered)
        else:
            list_msa_pb = list(x['sequence'] for x in list_msa_ordered)[:nummaxhits]

        return list_msa_pb
    
    def msa_split(self, list_msa_pb):
        """Split the processed MSA and return the list of commands for PSIBlast"""

        #Split for PSI blast
        if len(self.config_file.get('pssm', 'region_size')) > 0:
            #If we have a given region size then use that
            try:
                region_size = int(self.config_file.get('pssm', 'region_size'))
            except ValueError:
                print("[Homology Error] An incorrect value in the input for the region size, will be set as 20.")
                region_size = 20
        else:
            #Otherwise import a manual region 
            try:
                list_regions = literal_eval(self.config_file.get('pssm', 'manual_regions'))

                #Check to see if it is empty
                if len(list_regions) != 0:
                    list_regions = literal_eval(self.config_file.get('pssm', 'manual_regions'))
                    region_size = 0
                else:
                    print("[Homology Error] An incorrect value in the input for the manual regions, will be set as 20.")
                    region_size = 20
            except SyntaxError:
                print("[Homology Error] An incorrect value in the input for the manual regions, will be set as 20.")
                region_size = 20

        #Check our input
        if len(list_msa_pb) == 0:
            print("[msa_split Error] The input sequence list to be split is empty.")
            quit()

        #Create the regions or use a list
        if region_size == 0:
            #Use the input list
            regions = list_regions
        else:
            #Create the regions for a given size
            whole_regions = divmod(self.wtlen, region_size)[0]
            remainder = divmod(self.wtlen, region_size)[1]

            #Add one to the total number of regions if there is a remainder
            if remainder > 0:
                NumRegions = whole_regions + 1
            else:
                NumRegions = whole_regions

            #Calculate the low and high bounds of the regions
            regions = []
            for x in range(0, whole_regions + 1):

                if x < whole_regions:
                    low = x * region_size
                    high = (x * region_size) + region_size - 1
                else:
                    low = x * region_size
                    high = self.wtlen - 1

                #Create the list from the regions
                regions.append([low, high])

        #Inform the user
        print("[msa_split] Calculated Regions:")
        print(regions)
        print("[msa_split] Total number of regions: " + str(region_size)) 

        #Lets open the MSA (one sequence per line)
        list_psiblast_cmds = []
        for region in regions:
            start = region[0]
            end = region[1]
            counter = 0

            #Write the subalignment
            with open(self.directory + 'subalignment_' + self.proteinname + '_AAstart_' + str(start) + '.fasta', 'w') as outfile:
                for line in list_msa_pb:
                    subline = "".join(line[i].rstrip('\n') for i in range(start, end + 1))

                    #Remove sequences with insertions
                    if subline.find('-') == -1:
                        outfile.write(">" + str(counter) + "\n")
                        outfile.write(subline + "\n")
                        counter = counter + 1

            #Make a list of commands to run
            list_psiblast_cmds.append([
                "-subject",
                self.directory + self.proteinname + "_PSSM_WT.fasta",
                "-in_msa",
                self.directory + 'subalignment_' + self.proteinname + '_AAstart_' + str(start) + '.fasta',
                "-out_ascii_pssm",
                self.directory + self.proteinname + "_PSSM_" + str(start) + ".pssm",
                ])

        #Write a wildtype file for psiblast
        with open(self.directory + self.proteinname + "_PSSM_WT.fasta", "w") as file_wt:
            file_wt.write(">" + self.proteinname + "\n" + self.wtaa)

        return list_psiblast_cmds
    
    def pssm_file_import(self):
        """Combine the PSSM outputs into a lookup dict."""

        #Import our listdir function
        from os import listdir
    
        #Based on the PSSM output file. First two columns are wt and location.
        #The rest is the other columns with 0 = pssm score, 1 = percents
        pssm_table = [['NoneZ', None], ['NoneZ', None], 
        ['A', 0], ['R', 0], ['N', 0], ['D', 0], ['C', 0], ['Q', 0], ['E', 0], ['G', 0], ['H', 0], ['I', 0], 
        ['L', 0], ['K', 0], ['M', 0], ['F', 0], ['P', 0], ['S', 0], ['T', 0], ['W', 0], ['Y', 0], ['V', 0], 
        ['A', 1], ['R', 1], ['N', 1], ['D', 1], ['C', 1], ['Q', 1], ['E', 1], ['G', 1], ['H', 1], ['I', 1], 
        ['L', 1], ['K', 1], ['M', 1], ['F', 1], ['P', 1], ['S', 1], ['T', 1], ['W', 1], ['Y', 1], ['V', 1]]
        
        #Mutations matrix
        #[ResID][MutID[1]][0 = PSSMMetric, 1 = Percents]
        dict_pssm = {j: {i: [None, None] for i in self.aa_table} for j in range(1, self.wtlen + 1)}

        #Read all of the files
        for filename in listdir(self.directory):
            if filename.endswith(".pssm"):
                
                #Get the starting index number from the filename
                pssm_index = int(filename.split('_')[-1].split('.')[0])

                #Open the file again for reading
                with open(self.directory + filename, 'r') as file_pssm:
                    list_pssm_lines = file_pssm.readlines()

                #Parse each line
                for i in range(3, len(list_pssm_lines) - 4):

                    #Compact the line by removing non white space
                    list_line = list(c for c in list_pssm_lines[i].split(' ') if c != "")

                    #Add to our lookup table
                    for j in range(2, 42):
                        dict_pssm[int(list_line[0]) + pssm_index][pssm_table[j][0]][pssm_table[j][1]] = list_line[j]

        return dict_pssm

    def msa_freq(self, list_msa):
        """Process the msa and get the sitewise frequencies."""

        #Create a dictionary of site number : {'AA': 0, ...}
        dict_count = {i: {a:0 for a in self.aa_table} 
                      for i in range(1, len(self.wtaa) + 1)}

        #Iterate the msa list and add it to the dict_count
        for sequence in list_msa:
            for location in range(1, len(sequence) + 1):

                #Calculate the amino acid
                aa = sequence[location - 1]

                #Pass on the gaps '-'
                if aa not in self.aa_table:
                    continue
                
                #Add one
                dict_count[location][aa] = dict_count[location][aa] + 1

        #Get the residue totals
        dict_total = {i: sum([dict_count[i][a] for a in self.aa_table]) 
                      for i in range(1, len(self.wtaa) + 1)}

        #Get the residue frequencies
        dict_freq = {i: {a: dict_count[i][a]/dict_total[i] for a in self.aa_table} 
                      for i in range(1, len(self.wtaa) + 1)}

        return dict_freq

    def freq_output_heat(self, dict_freq):
        """This makes a CSV style heatmap."""
        
        #Print off the Number
        numbering = ' ,' + ','.join(str(i + 1) for i in range(0, self.wtlen))

        #Print off the WT Residue
        wtresi = ' ,' + ','.join(self.wtaa[i] for i in range(0, self.wtlen))

        #Setup the two new strings
        string_freq = ""

        #Loop each heatmap amino acid
        for amino in self.aa_table:

            #Add the amino acid to the start of the list
            string_freq = string_freq + amino + ","

            #Go through the len of the protein and construct the csv
            string_freq = string_freq + ','.join(str(dict_freq[j + 1][amino]) 
                                                 for j in range(0, self.wtlen)) + "\n"

        #Write the heatmap to a newfile
        with open(self.directory + self.proteinname + '_freq_Heatmap.csv', 'w') as outfile:
            outfile.write('\n'.join([
                numbering,
                wtresi,
                "Site Usage Frequencies",
                string_freq]))

        return "[output] Output saved as " + self.directory + self.proteinname + '_freq_Heatmap.csv'

    def pssm_output_heat(self, dict_pssm):
        """This makes a CSV style heatmap."""
        
        #Print off the Number
        numbering = ' ,' + ','.join(str(i + 1) for i in range(0, self.wtlen))

        #Print off the WT Residue
        wtresi = ' ,' + ','.join(self.wtaa[i] for i in range(0, self.wtlen))

        #Setup the two new strings
        string_pssm = ""
        string_percent = ""

        #Loop each heatmap amino acid
        for amino in self.aa_table:

            #Add the amino acid to the start of the list
            string_pssm = string_pssm + amino + ","
            string_percent = string_percent + amino + ","

            #Go through the len of the protein and construct the csv
            string_pssm = string_pssm + ','.join(str(dict_pssm[j + 1][amino][0]) 
                                                 for j in range(0, self.wtlen)) + "\n"
            
            string_percent = string_percent + ','.join(str(dict_pssm[j + 1][amino][1]) 
                                                       for j in range(0, self.wtlen)) + "\n"

        #Write the heatmap to a newfile
        with open(self.directory + self.proteinname + '_PSSM_Heatmap.csv', 'w') as outfile:
            outfile.write('\n'.join([
                numbering,
                wtresi,
                "Last position-specific scoring matrix computed",
                string_pssm,
                "Weighted observed percentages rounded down",
                string_percent]))

        return "[output] Output saved as " + self.directory + self.proteinname + '_PSSM_Heatmap.csv'

    def pssm_output_csv(self, dict_pssm):
        """Output our pssm dict as a csv column format."""
        
        #Write the dict to a newfile
        with open(self.directory + self.proteinname + '_PSSM_Columns.csv', 'w') as outfile:
            outfile.write("Location,Mutation,PSSM,Percent\n")
            outfile.write("\n".join([
                            ",".join([
                             str(location), 
                             amino, 
                             dict_pssm[location][amino][0], 
                             dict_pssm[location][amino][1]
                             ])
                            for location in dict_pssm 
                            for amino in dict_pssm[location]
                            ]))

        return "[output] Output saved as " + self.directory + self.proteinname + '_PSSM_Columns.csv'
    
    def xy_scatter(self, x, y):
        """Plot the dist of the read counts"""

        #Import the graphing toolkit
        from matplotlib import use, rcParams

        #Set the matplotlib base settings (headless, font)
        #use('Agg')
        rcParams['font.sans-serif']='Arial'
        rcParams['font.size']='18'
        rcParams.update({'figure.autolayout': True}) #Prevent labels from getting cut off

        from matplotlib import pyplot as plt

        #Check our axis limits
        if min(x) < self.x_min:
            self.x_min = min(x)

        if min(y) < self.y_min:
            self.y_min = min(y)

        if max(x) > self.x_max:
            self.x_max = max(x)

        if max(y) > self.y_max:
            self.y_max = max(y)

        #Create the scatterplot
        plt.scatter(x, y, s=1, color='black')

        #Set the axis limits
        plt.xlim(self.x_min, self.x_max)
        plt.ylim(self.y_min, self.y_max)
        plt.gca().set_aspect('equal', adjustable='box')
        
        plt.draw()

        #Label the graph
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        #If we are measuring regression then print it on the graph
        if self.bool_regression:
            #Import Pearsons Correlation
            from scipy.stats.stats import pearsonr

            #Calculate r
            float_pearsonr = round(pearsonr(x, y)[0], 2)

            #Add the legend
            textstr = "r = " + str(float_pearsonr)

            #Calculate top left
            text_x = self.x_max + ((self.x_max - self.x_min) / 20)
            text_y = (self.y_max - self.y_min) / 2
            plt.text(text_x, text_y, textstr, fontsize=18, 
                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))

        #Save the figure as a png
        plt.savefig(self.directory + self.proteinname + "_scatter.ps")
        plt.savefig(self.directory + self.proteinname + "_scatter.png")
        plt.show()

        return x, y

    def analysis_site_fit_homology_plot(self, dict_homology, dict_dataset, dataset_type):
        """
        Make a scatter plot of number of benefical versus highest homology variant
        Y = Number of benefical (i.e. #2SD) / total observed
        X = Highest homology variant (0.1 for frequency, -inf to inf for PSSM)
        """
        try:
            self.dataset_y = self.config_file.get('analysis_sitefitness_homology', 'dataset_y')
            
            self.y_column = self.config_file.get('analysis_sitefitness_homology', 'y_column').lower()
            self.y_threshold = float(self.config_file.get('analysis_sitefitness_homology', 'y_threshold'))

            self.x_label = self.config_file.get('analysis_sitefitness_homology', 'x_axis_label')
            self.y_label = self.config_file.get('analysis_sitefitness_homology', 'y_axis_label')

            self.x_min = float(self.config_file.get('analysis_sitefitness_homology', 'x_axis_min'))
            self.x_max = float(self.config_file.get('analysis_sitefitness_homology', 'x_axis_max'))
            self.y_min = float(self.config_file.get('analysis_sitefitness_homology', 'y_axis_min'))
            self.y_max = float(self.config_file.get('analysis_sitefitness_homology', 'y_axis_max'))

            #Set the output prefix
            self.output_prefix = self.directory + self.proteinname

            #Set the regression option
            self.bool_regression = bool(self.config_file.get('analysis_sitefitness_homology', 'regression'))
        except NoOptionError:
            print("[Error] Missing config file elements.")
            quit()

        #Make a blank dict as we need to save order
        x = {}
        y = {}

        #Loop the locations
        for location in dict_dataset[self.dataset_y]:

            #Setup a temp list
            temp_y = []

            #Loop the muts
            for mut in dict_dataset[self.dataset_y][location]:

                #Skip WT (Will need to correct for different starting points)
                if mut == dict_dataset[self.dataset_y][location][mut]['wt_residue']:
                    continue

                #Skip stops
                if mut == "*":
                    continue

                #If we have an observed log2
                if dict_dataset[self.dataset_y][location][mut][self.y_column] != "NaN":
                    temp_y.append(float(dict_dataset[self.dataset_y][location][mut][self.y_column]))

            #Add our tally of number of above threshold mutations
            count_threshold = sum([1 for mut in temp_y if mut >= self.y_threshold])
            y[location] = count_threshold / len(temp_y)

        #Loop the locations
        for location in dict_homology:

            #Temp count of highest value
            max_value = 0

            #Loop the mutations
            for mutation in dict_homology[location]:

                #If we have a higher value
                if dataset_type == "freq":
                    if dict_homology[location][mutation] > max_value:
                        max_value = dict_homology[location][mutation]
                else:
                    if int(dict_homology[location][mutation][0]) > max_value:
                        max_value = int(dict_homology[location][mutation][0])

            #Save our new max
            x[location] = max_value

        #Parse the dicts to return to lists
        x_list = []
        y_list = []
        for location in x:
            x_list.append(x[location])
            y_list.append(y[location])

        #Write the dict to a newfile
        with open(self.directory + self.proteinname + '_analysis_sitefitness_homology.csv', 'w') as outfile:
            outfile.write("Homology,Fraction_Over_Threshold\n")
            outfile.write("\n".join([
                            ",".join([str(x_list[i]), str(y_list[i])])
                            for i in range(0, len(x_list))
                            ]))

        #Make our graph
        self.xy_scatter(x_list, y_list)
        return

    def analysis_site_fit_homology_classifier(self, dict_homology, dict_dataset, dataset_type):
        """
        Make a scatter plot of number of benefical versus highest homology variant
        Y = Number of benefical (i.e. #2SD) / total observed
        X = Highest homology variant (0.1 for frequency, -inf to inf for PSSM)
        """
        try:
            self.dataset_y = self.config_file.get('analysis_sitefitness_homology', 'dataset_y')
            
            self.y_column = self.config_file.get('analysis_sitefitness_homology', 'y_column').lower()
            self.y_threshold = float(self.config_file.get('analysis_sitefitness_homology', 'y_threshold'))

            #Set the output prefix
            self.output_prefix = self.directory + self.proteinname
        except NoOptionError:
            print("[Error] Missing config file elements.")
            quit()

        #Make a blank dict as we need to save order
        x = {}
        y = {}

        #Loop the locations
        for location in dict_homology:

            #Temp count of highest value
            max_value = 0

            #Loop the mutations
            for mutation in dict_homology[location]:

                #If we have a higher value
                if dataset_type == "freq":
                    if dict_homology[location][mutation] > max_value:
                        max_value = dict_homology[location][mutation]
                else:
                    if int(dict_homology[location][mutation][0]) > max_value:
                        max_value = int(dict_homology[location][mutation][0])

            #Save our new max
            x[location] = max_value

        #Loop the locations
        for location in dict_dataset[self.dataset_y]:

            #Setup our dict
            y[location] = {"psd":0, "neutral":0, "msd":0}

            #Loop the muts
            for mut in dict_dataset[self.dataset_y][location]:

                #Skip WT (Will need to correct for different starting points)
                if mut == dict_dataset[self.dataset_y][location][mut]['wt_residue']:
                    continue

                #Skip stops
                if mut == "*":
                    continue

                #If we have an observed log2
                if dict_dataset[self.dataset_y][location][mut][self.y_column] == "NaN":
                    continue

                #Add to our counters
                if float(dict_dataset[self.dataset_y][location][mut][self.y_column]) >= self.y_threshold:
                    y[location]["psd"] = y[location]["psd"] + 1
                elif float(dict_dataset[self.dataset_y][location][mut][self.y_column]) <= (-1 * self.y_threshold):
                    y[location]["msd"] = y[location]["msd"] + 1
                else:
                    y[location]["neutral"] = y[location]["neutral"] + 1
                
        #Count and report
        dict_counts = {
            "total_1":0, "total_gp9":0, "total_gp75":0, "total_gp50":0,
            "total_gp25":0, "total_lp50":0, "total_lp25":0,

            "del_1":0, "del_gp9":0, "del_gp75":0, "del_gp50":0,
            "del_gp25":0, "del_lp50":0, "del_lp25":0,

            "neu_1":0, "neu_gp9":0, "neu_gp75":0, "neu_gp50":0,
            "neu_gp25":0, "neu_lp50":0, "neu_lp25":0,

            "ben_1":0, "ben_gp9":0, "ben_gp75":0, "ben_gp50":0,
            "ben_gp25":0, "ben_lp50":0, "ben_lp25":0,
            }
        
        #Loop the x dict
        for location in x:        

            if x[location] == 1:
                dict_counts["total_1"] = dict_counts["total_1"] + 1
                dict_counts["del_1"] = dict_counts["del_1"] + y[location]["msd"]
                dict_counts["neu_1"] = dict_counts["neu_1"] + y[location]["neutral"]
                dict_counts["ben_1"] = dict_counts["ben_1"] + y[location]["psd"]

            if x[location] >= 0.9:
                dict_counts["total_gp9"] = dict_counts["total_gp9"] + 1
                dict_counts["del_gp9"] = dict_counts["del_gp9"] + y[location]["msd"]
                dict_counts["neu_gp9"] = dict_counts["neu_gp9"] + y[location]["neutral"]
                dict_counts["ben_gp9"] = dict_counts["ben_gp9"] + y[location]["psd"]

            if x[location] >= 0.75:
                dict_counts["total_gp75"] = dict_counts["total_gp75"] + 1
                dict_counts["del_gp75"] = dict_counts["del_gp75"] + y[location]["msd"]
                dict_counts["neu_gp75"] = dict_counts["neu_gp75"] + y[location]["neutral"]
                dict_counts["ben_gp75"] = dict_counts["ben_gp75"] + y[location]["psd"]

            if x[location] >= 0.5:
                dict_counts["total_gp50"] = dict_counts["total_gp50"] + 1
                dict_counts["del_gp50"] = dict_counts["del_gp50"] + y[location]["msd"]
                dict_counts["neu_gp50"] = dict_counts["neu_gp50"] + y[location]["neutral"]
                dict_counts["ben_gp50"] = dict_counts["ben_gp50"] + y[location]["psd"]

            if x[location] >= 0.25:
                dict_counts["total_gp25"] = dict_counts["total_gp25"] + 1
                dict_counts["del_gp25"] = dict_counts["del_gp25"] + y[location]["msd"]
                dict_counts["neu_gp25"] = dict_counts["neu_gp25"] + y[location]["neutral"]
                dict_counts["ben_gp25"] = dict_counts["ben_gp25"] + y[location]["psd"]

            if x[location] <= 0.5:
                dict_counts["total_lp50"] = dict_counts["total_lp50"] + 1
                dict_counts["del_lp50"] = dict_counts["del_lp50"] + y[location]["msd"]
                dict_counts["neu_lp50"] = dict_counts["neu_lp50"] + y[location]["neutral"]
                dict_counts["ben_lp50"] = dict_counts["ben_lp50"] + y[location]["psd"]

            if x[location] <= 0.25:
                dict_counts["total_lp25"] = dict_counts["total_lp25"] + 1
                dict_counts["del_lp25"] = dict_counts["del_lp25"] + y[location]["msd"]
                dict_counts["neu_lp25"] = dict_counts["neu_lp25"] + y[location]["neutral"]
                dict_counts["ben_lp25"] = dict_counts["ben_lp25"] + y[location]["psd"]

        print(','.join(map(str, [
                "", "Fraction", "Total", "Del", "Neu", "Ben", "\n",                
                "1.0", dict_counts["total_1"],
                dict_counts["del_1"], dict_counts["neu_1"], dict_counts["ben_1"],
                "\n",
                ">= 0.9", dict_counts["total_gp9"],
                dict_counts["del_gp9"], dict_counts["neu_gp9"], dict_counts["ben_gp9"],
                "\n",
                ">= 0.75", dict_counts["total_gp75"],
                dict_counts["del_gp75"], dict_counts["neu_gp75"], dict_counts["ben_gp75"],
                "\n",
                ">= 0.5", dict_counts["total_gp50"],
                dict_counts["del_gp50"], dict_counts["neu_gp50"], dict_counts["ben_gp50"],
                "\n",
                ">= 0.25", dict_counts["total_gp25"],
                dict_counts["del_gp25"], dict_counts["neu_gp25"], dict_counts["ben_gp25"],
                "\n",
                "<= 0.5", dict_counts["total_lp50"],
                dict_counts["del_lp50"], dict_counts["neu_lp50"], dict_counts["ben_lp50"],
                "\n",
                "<= 0.25", dict_counts["total_lp25"],
                dict_counts["del_lp25"], dict_counts["neu_lp25"], dict_counts["ben_lp25"],
                "\n"])))



        #TODO - Make this for PSSMs

        return

    def classified_count_pssm_csv(self, dict_pssm, dict_classified, dict_merged_datasets, dataset_name, column):
        """Output a csv file of the fitness values"""

        """
        PSSM
        >= 3
        < 3 & >= 0
        < 0
        """

        #Get the unique class types by converting to a set
        set_class_types = set([dict_classified[loc][mut] 
                               for loc in dict_classified 
                               for mut in dict_classified[loc]])
        
        #Loop each class type
        for class_type in set_class_types:

            #Setup our lists
            list_threeplus = []
            list_zerotothree = []
            list_zerobelow = []

            #Loop the locations
            for loc in dict_classified:

                #Loop the mutations
                for mut in dict_classified[loc]:

                    #Only record from our class type
                    if dict_classified[loc][mut] != class_type:
                        continue

                    #Skip stops
                    if mut == "*":
                        continue

                    #>= 3
                    if int(dict_pssm[loc][mut][0]) >= 3:
                        list_threeplus.append(dict_merged_datasets[dataset_name][loc][mut][column])

                    #< 3 & >= 0
                    if int(dict_pssm[loc][mut][0]) >= 0 and int(dict_pssm[loc][mut][0]) < 3:
                        list_zerotothree.append(dict_merged_datasets[dataset_name][loc][mut][column])

                    #< 0
                    if int(dict_pssm[loc][mut][0]) < 0:
                        list_zerobelow.append(dict_merged_datasets[dataset_name][loc][mut][column])

            #Get our list lengths
            lena = len(list_threeplus)
            lenb = len(list_zerotothree)
            lenc = len(list_zerobelow)

            #Open the output file and write it
            with open(self.directory + self.proteinname + "_" 
                      + class_type + "_" 
                      + column + "_pssm_columns.csv", "w") as csv_out:

                #Write the header
                csv_out.write("PSSM_>=3,PSSM_>=0_<3,PSSM_<0\n")

                #Get the max length
                maxlen = max([lena, lenb, lenc])

                #Loop and write
                for i in range(0, maxlen + 1):
                
                    if i >= lena:
                        stra = ""
                    else:
                        stra = list_threeplus[i]

                    if i >= lenb:
                        strb = ""
                    else:
                        strb = list_zerotothree[i]

                    if i >= lenc:
                        strc = ""
                    else:
                        strc = list_zerobelow[i]

                    csv_out.write(str(stra) + "," + str(strb) + "," + str(strc) + "\n")

        return "Wrote CSV"

    def classified_count_pssm(self, dict_pssm, dict_classified):
        """Count our classifiers"""

        #Import our counter
        from collections import Counter

        """
        PSSM
        >= 3
        < 3 & >= 0
        < 0
        """

        #Setup our lists
        list_threeplus = []
        list_zerotothree = []
        list_zerobelow = []

        #Loop the locations
        for loc in dict_classified:

            #Loop the mutations
            for mut in dict_classified[loc]:

                #Skip stops
                if mut == "*":
                    continue

                #>= 3
                if int(dict_pssm[loc][mut][0]) >= 3:
                    list_threeplus.append(dict_classified[loc][mut])

                #< 3 & >= 0
                if int(dict_pssm[loc][mut][0]) >= 0 and int(dict_pssm[loc][mut][0]) < 3:
                    list_zerotothree.append(dict_classified[loc][mut])

                #< 0
                if int(dict_pssm[loc][mut][0]) < 0:
                    list_zerobelow.append(dict_classified[loc][mut])

        #Count the lists
        str_return = '\n'.join(map(str, [
            "PSSM: >= 3",
            pretty_counter_dicts(dict(Counter(list_threeplus))),
            "",
            "PSSM: < 3 & >= 0",
            pretty_counter_dicts(dict(Counter(list_zerotothree))),
            "",
            "PSSM: < 0",
            pretty_counter_dicts(dict(Counter(list_zerobelow)))
            ]))

        print(str_return)
        return str_return

if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[Homology Classifier Error] This classifier needs to be ran within the context of PACT.")
