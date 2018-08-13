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

"""Classifier Color - Assign colors to classifiers"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Classifier Color Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from pact.pact_common import get_bool

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class classifier_color:
    """Classifier Color - Assign colors to classifiers"""

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
            self.dataset = self.config_file.get('classifier_color', 'dataset')
            self.classifier = self.config_file.get('classifier_color', 'classifier')
            self.classifier_key = self.config_file.get('classifier_color', 'classifer_key').lower()
        except NoOptionError:
            print("[Classifier Color Error] Missing config file elements.")

        return

    def pdb_parse(self, dict_merged_datasets, dict_classifier):
        """Parse the pdb and return a dict for the coloring"""
        
        #The classifier line will be classifier: pdb, file.pdb
        pdb_file = self.config_file.get('classifier_color', 'pdb_file')
        pdb_chain = self.config_file.get('classifier_color', 'pdb_chain').upper()

        #Create a dict to work into
        dict_custom_color = {}

        #Loop our dataset
        for loc in dict_merged_datasets[self.dataset]:

            #Add our location to the dict
            dict_custom_color[loc] = {}

            #Loop our mutations
            for mut in dict_merged_datasets[self.dataset][loc]:

                #What classifer key to work with
                if self.classifier_key == "frac_burial":

                    #Get the info from the config file
                    burial_color = self.config_file.get('classifier_color', 'burial_color').lower()
                    burial_value = float(self.config_file.get('classifier_color', 'burial_value'))
                    burial_equality = self.config_file.get('classifier_color', 'burial_equality')
                    burial_othercolor = self.config_file.get('classifier_color', 'burial_othercolor').lower()

                    #May have to include chain

                    #Verify that the location is present in our structure
                    if loc not in dict_classifier[pdb_file]['dssp'][pdb_chain]:                        
                        dict_custom_color[loc][mut] = "black"
                        continue

                    #Test and append to our dict
                    if get_bool(dict_classifier[pdb_file]['dssp'][pdb_chain][loc]['frac_burial'], burial_equality, burial_value):
                        dict_custom_color[loc][mut] = burial_color
                    else:
                        dict_custom_color[loc][mut] = burial_othercolor

        return dict_custom_color

    def classifier_color(self, dict_merged_datasets, dict_classifier, classifier_type):
        """Main entrypoint"""

        #Decide on type of classifier
        if classifier_type == "pdb":

            #Call the pdb parser
            dict_custom_color = self.pdb_parse(dict_merged_datasets, dict_classifier)

        return dict_custom_color

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Classifier Color Error] This file needs to be ran within the context of PACT.")
