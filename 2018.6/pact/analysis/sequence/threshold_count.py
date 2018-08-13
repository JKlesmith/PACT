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

"""Threshold Count - Count the number of mutations above or between sets of values"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Threshold Count Error] Your Python interpreter is too old."
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

class threshold_count:
    """Count the number of mutations above or between sets of values"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""

        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("threshold_count"):
            print("[Threshold Count Error] the config file is missing the section [threshold_count]")
            quit()
        
        #Classifier specific variables
        try:
            self.dataset = self.config_file.get('threshold_count', 'dataset')
            self.column = self.config_file.get('threshold_count', 'column').lower()
            self.cutoff = float(self.config_file.get('threshold_count', 'cutoff'))
        except NoOptionError:
            print("[Threshold Count Error] Missing config file elements.")
            quit()

        return

    def parse_dataset(self, dict_merged_datasets):
        """Parse the fitness dict"""

        #Check to see if our dataset is in our merged datasets
        if self.dataset not in dict_merged_datasets:
            print("[Threshold Count Error] Check the config file, the dataset wanted was not combined.")
            quit()
        
        #Make a empty list
        count_above_cut = 0
        count_below_cut = 0

        #Parse the dict_fitness
        for key in dict_merged_datasets[self.dataset]:

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset][key]:

                #If we have a fitness
                if dict_merged_datasets[self.dataset][key][mut][self.column] == "NaN":
                    continue

                #If we have a fitness
                if (dict_merged_datasets[self.dataset][key][mut]['wt_residue'] == mut):
                    continue
                
                #If our column is above our cutoff
                if dict_merged_datasets[self.dataset][key][mut][self.column] >= self.cutoff:
                    count_above_cut = count_above_cut + 1
                else:
                    count_below_cut = count_below_cut + 1

        return count_above_cut, count_below_cut

    def threshold_count(self, dict_merged_datasets):
        """Run the workflow to count the mutations above a threshold"""

        #Parse the dataset
        count_above_cut, count_below_cut = self.parse_dataset(dict_merged_datasets)

        output_string = '\r\n'.join(map(str, [
        "Threshold Value: ",
        self.cutoff,
        "Column: ",
        self.column,
        "Dataset: ",
        self.dataset,
        "Above Cutoff: ",
        count_above_cut,
        "Below Cutoff: ",
        count_below_cut]))
            
        print(output_string)

        return output_string

if __name__ == '__main__':

    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[Threshold Count Error] This classifier needs to be ran within the context of PACT.") 
