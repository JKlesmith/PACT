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

"""aa_fitmet_compare - compare two groups of amino acids"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[aa_fitmet_compare Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError, NoSectionError

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class aa_fitmet_compare:
    """Compare the fitness of two groups of amino acids"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""

        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("aa_compare_ttest"):
            print("[Amino Acid Compare Error] the config file is missing the section [aa_compare_ttest]")
            quit()

        #Classifier specific variables
        self.aa_table = '*FWYPMILVAGCSTNQDEHKR'
        self.directory = self.dict_protocolconfig['directory']
        self.output_prefix = self.config_file.get('global', 'output_prefix')

        #Specific Options
        try:
            self.dataset = self.config_file.get('aa_compare_ttest', 'dataset')
            
            self.group_a = self.config_file.get('aa_compare_ttest', 'group_a').upper()
            self.group_b = self.config_file.get('aa_compare_ttest', 'group_b').upper()
            
            self.group_a_title = self.config_file.get('aa_compare_ttest', 'group_a_title')
            self.group_b_title = self.config_file.get('aa_compare_ttest', 'group_b_title')

            self.y_axis_label = self.config_file.get('aa_compare_ttest', 'y_axis_label')

            self.column = self.config_file.get('aa_compare_ttest', 'column').lower()

            self.excludewt = self.config_file.get('aa_compare_ttest', 'exclude_wt').lower()

            self.pact_headless = self.config_file.get('aa_compare_ttest', 'headless').lower()

        except NoOptionError:
            print("[Amino Acid Compare Error] Missing config file elements.")
            quit()

        #Check the validity of the groups
        for letter in self.group_a:
            if letter not in self.aa_table:
                print("[Amino Acid Compare Error] Unknown letter in group a.")
                quit()

        for letter in self.group_b:
            if letter not in self.aa_table:
                print("[Amino Acid Compare Error] Unknown letter in group b.")
                quit()

        return

    def parse_dictfit(self, dict_merged_datasets):
        """Parse the fitness dict"""

        #Check to see if our dataset is in our merged datasets
        if self.dataset not in dict_merged_datasets:
            print("[Amino Acid Compare Error] Check the config file, the dataset wanted was not combined.")
            quit()
        
        #Make a empty list
        list_grpa = []
        list_grpb = []

        #Parse the dict_fitness
        for key in dict_merged_datasets[self.dataset]:

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset][key]:

                #If we have a fitness
                if dict_merged_datasets[self.dataset][key][mut][self.column] == "NaN":
                    continue

                #If we have a fitness
                if (dict_merged_datasets[self.dataset][key][mut]['wt_residue'] == mut and self.excludewt == "true"):
                    continue

                #Set our bools
                bool_grpa = False

                #Parse the mutations (*,A,C)
                for letter in dict_merged_datasets[self.dataset][key][mut]['mutation'].split(","):

                    #Assign a group (group A will trump group B)
                    if letter in self.group_a:
                        bool_grpa = True
                
                #What is our group?
                if bool_grpa:
                    list_grpa.append(dict_merged_datasets[self.dataset][key][mut][self.column])
                else:
                    list_grpb.append(dict_merged_datasets[self.dataset][key][mut][self.column])

        return list_grpa, list_grpb
    
    def aa_dist_ttest(self, list_grpa, list_grpb):
        """Run a t-test of the two distributions and return the statistics."""

        #Import the t-test from scipy
        from scipy.stats import ttest_ind

        #Perform the t-test on the two lists
        stats, pvalue = ttest_ind(list_grpa, list_grpb, equal_var = False)

        #Print our info
        print("p-value")
        print(pvalue)
        print("N Group A")
        print(len(list_grpa))
        print("N Group B")
        print(len(list_grpb))

        return stats, pvalue

    def create_boxplot(self, list_grpa, list_grpb, tvalue, pvalue):
        """Create a boxplot figure of the two distributions."""

        #Import the graphing toolkit
        from matplotlib import use, rcParams

        #Set the matplotlib base settings (headless, font)
        if self.pact_headless == "true":
            use('Agg')
        rcParams['font.sans-serif']='Arial'

        from matplotlib import pyplot as plt

        # Create a figure instance
        fig = plt.figure(1, figsize=(3, 6))

        # Create an axes instance
        ax = fig.add_subplot(111)

        # Create the boxplot
        bp = ax.boxplot([list_grpa, list_grpb], widths=(0.5, 0.5), sym='+')

        # Custom x-axis labels
        ax.set_xticklabels([self.group_a_title + "\n N=" + str(len(list_grpa)), 
                            self.group_b_title + "\n N=" + str(len(list_grpb))])

        #Add the legend
        textstr = "Welch's t-test\nt=" + str(tvalue) + "\np=" + str(pvalue)
        plt.text(1.1, 0.5, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')

        #Label the graph
        plt.ylabel(self.y_axis_label)
        
        #Save the figure as a png
        plt.savefig(self.directory + self.output_prefix + "_grp_vs_grp.png", dpi=600, bbox_inches='tight')

        if self.pact_headless != "true":
            plt.show()

        return

    def create_csv(self, list_grpa, list_grpb):
        """"Create a csv of the two groups"""

        lena = len(list_grpa)
        lenb = len(list_grpb)

        #Open the output file and write it
        with open(self.directory + self.output_prefix + "_grp_vs_grp.csv", "w") as csv_out:

            #Write the header
            csv_out.write(self.group_a_title + "," + self.group_b_title + "\n")

            #Get the max length
            if lena > lenb:
                maxlen = lena
            else:
                maxlen = lenb

            #Loop and write
            for i in range(0, maxlen + 1):
                
                if i >= lena:
                    stra = ""
                else:
                    stra = list_grpa[i]

                if i >= lenb:
                    strb = ""
                else:
                    strb = list_grpb[i]

                csv_out.write(str(stra) + "," + str(strb) + "\n")

        return

    def aa_fitmet_compare(self, dict_merged_datasets):
        """Run the workflow to compare groups of amino acids"""

        #Parse the fitness dict
        list_grpa, list_grpb = self.parse_dictfit(dict_merged_datasets)

        #Do a Welch's T-Test
        tvalue, pvalue = self.aa_dist_ttest(list_grpa, list_grpb)

        #Make the boxplot
        self.create_boxplot(list_grpa, list_grpb, tvalue, pvalue)

        #Output the CSV
        self.create_csv(list_grpa, list_grpb)

        #Print stats
        str_stats = '\n'.join([
            "Welch's T-Test:",
            "tvalue = " + str(tvalue),
            "pvalue = " + str(pvalue),
            "Wrote CSV group data"
            ])
        print(str_stats)

        return str_stats

if __name__ == '__main__':

    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[Amino Acid Compare Error] This classifier needs to be ran within the context of PACT.") 
