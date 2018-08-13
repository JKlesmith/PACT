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

"""dataset vs dataset reporter - analyze and plot two pact datasets against each other"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Set vs Set Reporter Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from pact.pact_common import file_checker
from scipy.stats import linregress
from statistics import mean, stdev

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class set_vs_set_counter:
    """dataset vs dataset counter - analyze and plot two pact datasets against each other"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Set the output prefix
        self.directory = self.dict_protocolconfig['directory']
        self.output_prefix = self.config_file.get('global', 'output_prefix')

        #Options
        try:
            self.dataset_x = self.config_file.get('setvsset', 'dataset_x')
            self.dataset_y = self.config_file.get('setvsset', 'dataset_y')
            
            self.class_column = self.config_file.get('variant_classification', 'class_column').lower()
            self.class_threshold = float(self.obj_cfgparser.get("variant_classification", "class_threshold"))

        except NoOptionError:
            print("[Set vs Set Reporter Error] Missing config file elements.")
            quit()

        return

    def xy_list(self, dict_merged_datasets):
        """Return two lists of x and y"""

        #Make a blank list
        x = []
        y = []
        x_nan = []
        y_nan = []

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    x_nan.append(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]))
                    y_nan.append(float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Add our mutation
                x.append(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]))
                y.append(float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))

        return x, y, x_nan, y_nan




    def create_csv(self, x, y):
        """Output a CSV for other programs"""

        #Open the output file and write it
        with open(self.directory + self.output_prefix + "_x_vs_y.csv", "w") as csv_out:

            #Write the header
            csv_out.write(self.x_label + "," + self.y_label + "\n")

            #Create the output lines
            csv_out.write('\n'.join([str(xlist) + "," + str(ylist) for xlist,ylist in zip(x, y)]))
        return







    def set_vs_set_reporter(self, dict_merged_datasets):
        """Main entrypoint"""

        #Check to see if our dataset is in our merged datasets
        if (self.dataset_x not in dict_merged_datasets or
            self.dataset_y not in dict_merged_datasets):
            print("[Set vs Set Reporter Error] Check the config file, the dataset wanted was not combined.")
            quit()




        #Output a CSV of the data
        self.create_csv(x, y)



        #Do we calculate shared winners counts
        if self.count_shared == "true":
            #Get our basal counts of shared winners
            countx, county, countxy = self.shared_mutations(dict_merged_datasets)

            #Print shared winner stats
            str_shared = '\n'.join([
                "Dataset vs Dataset",
                "Dataset X: " + self.dataset_x,
                "Dataset Y: " + self.dataset_y,
                "Counts of nonsynonymous mutations above threshold: " + str(self.winner_threshold),
                "Dataset X: " + str(countx),
                "Dataset Y: " + str(county),
                "Shared: " + str(countxy),
                ])
            print(str_shared)
        else:
            str_shared = ""

        return ""

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Set vs Set Counter Error] This file needs to be ran within the context of PACT.")
