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

"""convert_csv_to_pact- Import a CSV then convert to pact format"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Read CSV Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError, NoSectionError
from pact.pact_common import file_checker, save_pact_file

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class convert_csv_to_pact:
    """Import fitness from a csv file"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("convert_csv_to_pact"):
            print("[Convert CSV to PACT Error] the config file is missing the section [convert_csv_to_pact]")
            quit()

        #Set the number of datasets
        try:
            self.numdatasets = int(self.config_file.get('convert_csv_to_pact', 'numdatasets'))
        except ValueError:
            print("[Convert CSV to PACT Error] The number of datasets is not set properly in the config file.")
        except NoOptionError:
            print("[Convert CSV to PACT Error] The convert_csv_to_pact config section is incorrect.")
            print("[Convert CSV to PACT Error] Missing: numdatasets")
            quit()

        #Get the dataset names
        self.dict_combined = {}

        try:
            for i in range(1, self.numdatasets + 1):
                self.dict_combined[self.config_file.get('convert_csv_to_pact', 'dataset_' + str(i))] = {}
        except NoOptionError:
            print("[Read CSV Error] The convert_csv_to_pact config file is incorrect.")
            print("[Read CSV Error] There is something wrong with the name of a option flag.")
            quit()

        return

    def read_csv_fitness(self):
        """From the config file it loads and combines pact files"""
        
        #For each dataset parse the config file for the files
        for dataset in self.dict_combined:

            #Info the user
            print("[Read CSV] Merging " +  dataset)

            #Test if the section is present
            if not self.config_file.has_section(dataset):
                print("[Read CSV Error] the config file is missing the section for " + dataset)

            #Load all of the files into a dict
            try:
                dict_section = dict(self.config_file.items(dataset))
            except NoSectionError:
                print("[Read CSV] Can not find the dataset given, check config file.")
                quit()

            #See if the file exists
            if file_checker(dict_section['file']):

                #Load the pact file
                with open(dict_section['file'], 'r') as infile:
                    list_file = infile.readlines()

                #Figure out which column is which
                column_location = list_file[0].split(',').index(dict_section['location'])

                #Figure out which column is which
                column_mutation = list_file[0].split(',').index(dict_section['mutation'])

                #Figure out which column is which
                column_fitness = list_file[0].split(',').index(dict_section['fitness'])

                #Figure out which column is which
                starting_index = int(dict_section['starting_index'])

                if starting_index == 0:
                    wt_offset = 0
                elif starting_index == 1:
                    wt_offset = -1

                #Add to our dict, skipping the header
                for line in list_file[1:]:

                    #Set our variables
                    loc = int(line.split(',')[column_location])
                    mut = line.split(',')[column_mutation]

                    #Add our location if not existing
                    if loc not in self.dict_combined[dataset]:
                        self.dict_combined[dataset][loc] = {}

                    #Warn the user if the key exists
                    if mut in self.dict_combined[dataset][loc]:
                        print("[Read CSV Warning] Key: " + str(loc) + str(mut) + " exists in dataset " + dataset)

                    #Add to the global dict
                    fitness = line.split(',')[column_fitness]

                    if fitness == 'NS':
                        fitness = 'NaN'
                    else:
                        fitness = float(fitness)

                    #Get the wild-type residue
                    wt_resi = self.config_file.get(dataset, 'wtaa')[loc + wt_offset]

                    self.dict_combined[dataset][loc][mut] = {'location':loc, 
                                                             'fitness':fitness, 
                                                             'mutation':mut, 
                                                             'wt_residue':wt_resi,
                                                             'sd_from_wt':'NaN'}

            #Save a pact file for our dataset
            save_pact_file(self.dict_combined[dataset], self.dict_protocolconfig['directory'] + '/' + dataset)

        return self.dict_combined

if __name__ == '__main__':
    #Remind the user that this file needs to be ran within the context of PACT
    print("[Read CSV Error] This file needs to be ran within the context of PACT.")
