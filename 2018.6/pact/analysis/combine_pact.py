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

"""combine_pact.py - Combine datasets"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Combine PACT Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError, NoSectionError
from pact.pact_common import file_checker
from pickle import load, UnpicklingError

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class combine_pact:

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("combinepact"):
            print("[Combine PACT Error] the config file is missing the section [combinepact]")
            quit()

        #Set the number of datasets
        try:
            self.numdatasets = int(self.config_file.get('combinepact', 'numdatasets'))
        except ValueError:
            print("[Combine Pact Error] The number of datasets is not set properly in the config file.")
        except NoOptionError:
            print("[Combine Pact Error] The combinepact config section is incorrect.")
            print("[Combine Pact Error] Missing: numdatasets")
            quit()

        #Get the dataset names
        self.dict_combined = {}

        try:
            for i in range(1, self.numdatasets + 1):
                self.dict_combined[self.config_file.get('combinepact', 'dataset_' + str(i))] = {}
        except NoOptionError:
            print("[Combine Pact Error] The combinepact config file is incorrect.")
            print("[Combine Pact Error] There is something wrong with the name of a option flag.")
            quit()

        return

    def combine_pact(self):
        """From the config file it loads and combines pact files"""
        
        #For each dataset parse the config file for the files
        for dataset in self.dict_combined:

            #Info the user
            print("[Combine Pact] Merging " +  dataset)

            #Test if the section is present
            if not self.config_file.has_section(dataset):
                print("[Combine Pact Error] the config file is missing the section for " + dataset)

            #Load all of the files into a dict
            try:
                dict_files = dict(self.config_file.items(dataset))
            except NoSectionError:
                print("[Combine PACT] Can not find the dataset given, check config file.")
                quit()

            #Load the file
            for file in dict_files:

                #See if the file exists
                if file_checker(dict_files[file]):

                    #Load the pact file
                    with open(dict_files[file], 'rb') as infile:
                        try:
                            pact_dict = load(infile)
                        except UnpicklingError:
                            print("[PACT Error] Check your .pact filename.")
                            quit()

                    #Add the pact file to the dataset dict
                    for key in pact_dict:

                        #Warn the user if the key exists
                        if key in self.dict_combined[dataset]:
                            print("[Combine Pact Warning] Key: " + str(key) + " exists in dataset " + dataset)

                        #Add to the global dict
                        self.dict_combined[dataset][key] = pact_dict[key]

        return self.dict_combined

if __name__ == '__main__':
    #Remind the user that this file needs to be ran within the context of PACT
    print("[Combine Pact Error] This file needs to be ran within the context of PACT.")
