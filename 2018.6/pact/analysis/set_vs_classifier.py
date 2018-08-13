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

"""dataset vs classifier"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Set vs Set Error] Your Python interpreter is too old. Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from pact.pact_common import file_checker, get_bool

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class set_vs_classifier:

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Options
        try:
            self.dataset_y = self.config_file.get('setvsclassifier', 'dataset_y')
            
            self.y_column = self.config_file.get('setvsclassifier', 'y_column').lower()

            self.x_label = self.config_file.get('setvsclassifier', 'x_axis_label')
            self.y_label = self.config_file.get('setvsclassifier', 'y_axis_label')

            self.x_min = float(self.config_file.get('setvsclassifier', 'x_axis_min'))
            self.x_max = float(self.config_file.get('setvsclassifier', 'x_axis_max'))
            self.y_min = float(self.config_file.get('setvsclassifier', 'y_axis_min'))
            self.y_max = float(self.config_file.get('setvsclassifier', 'y_axis_max'))

            #Set the output prefix
            self.output_prefix = self.config_file.get('global', 'output_prefix')

        except NoOptionError:
            print("[Set vs Set Error] Missing config file elements.")
            quit()

        return

    def y_list(self):
        """Return two lists of x and y"""

        #Make a blank list
        y = []
        color = []
        y_nan = []

        #Loop the locations
        #try:
        for key in self.dict_merged[self.dataset_y]:

            #Loop the muts
            for mut in self.dict_merged[self.dataset_y][key]:
                
                #Skip WT (Will need to correct for different starting points)
                if mut == self.dict_merged[self.dataset_y][key][mut]['wt_residue']:
                    continue

                if (self.dict_merged[self.dataset_y][key][mut][self.y_column] != "NaN"):
                    if (self.dict_merged[self.dataset_y][key][mut][self.y_column] > 2):
                        color.append(self.dict_classifier[key]['frac_burial'])
                    #y_nan.append(float(self.dict_merged[self.dataset_y][key][mut][self.y_column]))
                #else:
                    #y.append(float(self.dict_merged[self.dataset_y][key][mut][self.y_column]))

                    #Make a color list
                    #color.append(self.dict_classifier[key]['frac_burial'])
        #except KeyError:
        #    print("[Set vs Set Error] Check your config file as the x dataset name is not found.")
        #    quit()

        return y, color, y_nan

    def sety_hist(self):
        """Plot the dist of the read counts"""

        #Import the graphing toolkit
        from matplotlib import use, rcParams

        #Set the matplotlib base settings (headless, font)
        #use('Agg')
        #rcParams['font.sans-serif']='Arial'

        from matplotlib import pyplot as plt
        from numpy import ones_like, linspace

        #Get the lists of values
        y, color, y_nan = self.y_list()

        #Make the x-axis bins and y-axis weights (to make it out of 1)
        bins = linspace(0, 1, 100)
        weights = ones_like(color)/float(len(color)) #Weight the y-axis values to make it out of 1

        #Create the histogram
        fig, ax = plt.subplots()
        counts, bins, patches = ax.hist(color, bins, weights=weights, color='black')

        #Label the graph
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
    
        #Add the legend
        #textstr = "Mean: "+str(count_mean)+"\nMedian: "+str(count_median)+"\nN: "+str(len(list_counts))
        #plt.text(0.8, 0.95, textstr, transform=ax.transAxes, fontsize=12, 
        #         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))

        #Save the figure as a png
        plt.savefig(self.config_file.get('global', 'directory') + self.output_prefix + "_hist.png")
        plt.show()

        return y

    def xy_list(self):
        """Return two lists of x and y"""

        """Like to make a xy scatter of fitness or sd vs frac buried (show that best ones is from surface)"""

        #Make a blank list
        x = []
        y = []
        color = []


        #Loop the locations
        #try:
        for key in self.dict_merged[self.dataset_y]:

            #Loop the muts
            for mut in self.dict_merged[self.dataset_y][key]:
                
                #Skip WT (Will need to correct for different starting points)
                if mut == self.dict_merged[self.dataset_y][key][mut]['wt_residue']:
                    continue

                if (self.dict_merged[self.dataset_y][key][mut][self.y_column] != "NaN"):
                    #if (self.dict_merged[self.dataset_y][key][mut][self.y_column] > 2):

                    x.append(self.dict_classifier[key]['frac_burial'])
                    y.append(float(self.dict_merged[self.dataset_y][key][mut][self.y_column]))
                    color.append(self.dict_classifier[key]['classifer_color'])

        #except KeyError:
        #    print("[Set vs Set Error] Check your config file as the x dataset name is not found.")
        #    quit()

        return x, y, color

    def xy_scatter(self):
        """Plot the dist of the read counts"""

        #Import the graphing toolkit
        from matplotlib import use, rcParams

        #Set the matplotlib base settings (headless, font)
        #use('Agg')
        #rcParams['font.sans-serif']='Arial'

        from matplotlib import pyplot as plt

        #Get the two lists of values
        x, y, color = self.xy_list()

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
        plt.scatter(x, y, s=1, color=color)

        #Set the axis limits
        plt.xlim(self.x_min, self.x_max)
        plt.ylim(self.y_min, self.y_max)
        #plt.gca().set_aspect('equal', adjustable='box')
        #plt.draw()

        #Label the graph
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        #Plot the origin lines
        plt.axhline(0, color='black', linestyle='--')
        plt.axvline(0, color='black', linestyle='--')


        #Save the figure as a png
        plt.savefig(self.config_file.get('global', 'directory') + self.output_prefix + "_scatter.png")
        plt.show()

        return x, y

    def chi_square_frac_burial(self, frac_buried, relation, std_dev):
        """Quantify the number of mutation changes that are bene, neutral, or del at a given fraction buried"""

        print("[Analysis] Fraction Buried: " + relation + str(frac_buried))

        #Setup our couters
        count_total_deleterious = 0
        count_total_neutral = 0
        count_total_beneficial = 0

        count_deleterious_hptopc = 0
        count_neutral_hptopc = 0
        count_beneficial_hptopc = 0
        
        count_deleterious_pctohp = 0
        count_neutral_pctohp = 0
        count_beneficial_pctohp = 0 
        
        count_deleterious_hptohp = 0
        count_neutral_hptohp = 0
        count_beneficial_hptohp = 0

        count_deleterious_pctopc = 0
        count_neutral_pctopc = 0
        count_beneficial_pctopc = 0

        count_deleterious_xtopro = 0
        count_neutral_xtopro = 0
        count_beneficial_xtopro = 0

        count_deleterious_xtocys = 0
        count_neutral_xtocys = 0
        count_beneficial_xtocys = 0 

        #Loop the residues
        for key in self.dict_merged[self.dataset_y]:

            #Loop the muts
            for mut in self.dict_merged[self.dataset_y][key]:
                
                #Skip WT
                if mut == self.dict_merged[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip not a number
                if self.dict_merged[self.dataset_y][key][mut][self.y_column] == "NaN":
                    continue

                #Skip stop codons
                if mut == "*":
                    continue

                #Setup our table
                if self.get_bool(self.dict_structure_classifier[key]['frac_burial'], relation, frac_buried):

                    #Setup a >=2 polar/charge to HP list
                    if self.dict_merged[self.dataset_y][key][mut][self.y_column] >= std_dev:
                        count_total_beneficial = count_total_beneficial + 1

                        if (self.dict_sequence_classifier[key][mut]['classifer'] == "LB" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PB" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CB"):

                            count_beneficial_pctohp = count_beneficial_pctohp + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BL" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PL" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CL"):
                        
                            count_beneficial_hptopc = count_beneficial_hptopc + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "LL"):
                            count_beneficial_pctopc = count_beneficial_pctopc + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BB"):
                            count_beneficial_hptohp = count_beneficial_hptohp + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BP" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "LP" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CP"):
                        
                            count_beneficial_xtopro = count_beneficial_xtopro + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BC" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "LC" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PC"):
                        
                            count_beneficial_xtocys = count_beneficial_xtocys + 1

                    elif self.dict_merged[self.dataset_y][key][mut][self.y_column] <= -1 * std_dev:
                        count_total_deleterious = count_total_deleterious + 1

                        if (self.dict_sequence_classifier[key][mut]['classifer'] == "LB" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PB" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CB"):

                            count_deleterious_pctohp = count_deleterious_pctohp + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BL" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PL" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CL"):
                        
                            count_deleterious_hptopc = count_deleterious_hptopc + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "LL"):
                            count_deleterious_pctopc = count_deleterious_pctopc + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BB"):
                            count_deleterious_hptohp = count_deleterious_hptohp + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BP" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "LP" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CP"):
                        
                            count_deleterious_xtopro = count_deleterious_xtopro + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BC" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "LC" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PC"):
                        
                            count_deleterious_xtocys = count_deleterious_xtocys + 1

                    else:
                        #Setup the neutral counts
                        count_total_neutral = count_total_neutral + 1

                        if (self.dict_sequence_classifier[key][mut]['classifer'] == "LB" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PB" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CB"):

                            count_neutral_pctohp = count_neutral_pctohp + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BL" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PL" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CL"):
                        
                            count_neutral_hptopc = count_neutral_hptopc + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "LL"):
                            count_neutral_pctopc = count_neutral_pctopc + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BB"):
                            count_neutral_hptohp = count_neutral_hptohp + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BP" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "LP" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "CP"):
                        
                            count_neutral_xtopro = count_neutral_xtopro + 1

                        elif (self.dict_sequence_classifier[key][mut]['classifer'] == "BC" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "LC" or
                            self.dict_sequence_classifier[key][mut]['classifer'] == "PC"):
                        
                            count_neutral_xtocys = count_neutral_xtocys + 1

        print('\t'.join(map(str, [
        "\tTotal Counts: Del Neu Ben\n",
        count_total_deleterious,
        count_total_neutral,
        count_total_beneficial,
        "\n",
        "Non-Hydrophobic to Hydrophobic: Del Neu Ben\n",
        count_deleterious_pctohp,
        count_neutral_pctohp,
        count_beneficial_pctohp, 
        "\n",
        "Non-Polar/Charged to Polar/Charged: Del Neu Ben\n",
        count_deleterious_hptopc,
        count_neutral_hptopc,
        count_beneficial_hptopc,
        "\n",
        "Polar/Charged to Polar/Charged: Del Neu Ben\n",
        count_deleterious_pctopc,
        count_neutral_pctopc,
        count_beneficial_pctopc, 
        "\n",
        "Hydrophobic to Hydrophobic: Del Neu Ben\n",
        count_deleterious_hptohp,
        count_neutral_hptohp,
        count_beneficial_hptohp,
        "\n",
        "X to Pro: Del Neu Ben\n",
        count_deleterious_xtopro,
        count_neutral_xtopro,
        count_beneficial_xtopro,
        "\n",
        "X to Cys: Del Neu Ben\n",
        count_deleterious_xtocys,
        count_neutral_xtocys,
        count_beneficial_xtocys,
        "\n"])))
        return

    def set_vs_classifier(self, dict_merged_datasets):
        """What the class does"""

        #Perform a regular hist
        #y = self.sety_hist()
        #self.xy_scatter()
        self.chi_square_frac_burial(0.0, '>=', 2)
        self.chi_square_frac_burial(0.9, '>=', 2)
        self.chi_square_frac_burial(0.85, '>=',2)
        self.chi_square_frac_burial(0.5, '<=', 2)

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Set vs Classifier Error] This file needs to be ran within the context of PACT.")

"""        
*****************************************
Set vs Classifier Section
*****************************************
            
if self.dict_workflow['setvsclassifier']:

    #Check to see if the section is there
    if not self.obj_cfgparser.has_section('setvsclassifier'):           
        print("[Protocols:Analysis Error] The setvsclassifier config file is incorrect.")
        print("[Protocols:Analysis Error] There is something wrong with the [setvsclassifier] section.")
        quit()

    #Import our setvsclassifier class
    try:
        from pact.analysis.set_vs_classifier import set_vs_classifier
    except ImportError:
        print("[Protocols:Analysis Error] set_vs_set was not found.")

    #Create the object then call the merger
    obj_svc = set_vs_classifier(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

    #Run the main routine
    obj_svc.set_vs_classifier(dict_merged_datasets)
"""