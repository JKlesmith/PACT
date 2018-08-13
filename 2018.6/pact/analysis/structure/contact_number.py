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

"""contact_number - Calculate the number of Ca/b atoms in a sphere around a Ca/b atom"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Contact Number Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from math import sqrt, pow
from pact.pact_common import pretty_counter_dicts, euc_dist

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class contact_number:
    """Class with methods to calculate the number of Ca or Cb atoms in a sphere around a Ca or Cb"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Setup the globals"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Get any config file entries
        if self.config_file.has_section("contact_number"):
            try:
                #Classifier specific variables
                self.pdb_file = self.config_file.get('contact_number', 'pdb_file')
                self.config_atoms = self.config_file.get('contact_number', 'atoms')
                self.config_chains = self.config_file.get('contact_number', 'chains')
                self.report_chain = self.config_file.get('contact_number', 'report_chain')
                self.config_distance = self.config_file.get('contact_number', 'distance')
            except NoOptionError:
                print("[Contact Number Error] the config file is missing any of the options: pdb_file, atoms, chains, or distance.")
                self.config_atoms = "CA"
                self.config_chains = "A"
                self.config_distance = 10
                self.report_chain = "A"
        else:
            print("[Contact Number Error] the config file is missing the section [contact_number]")

        return

    def contact_number(self, dict_pdb, key_pdb='', atom_include='', distance='', chain_include=''):
        """Calculate the contact number of a residue"""

        #Use the config file or the function call
        if len(atom_include) > 0:
            list_atoms = atom_include.upper().split(',')
        else:
            list_atoms = self.config_atoms.upper().split(',')

        if len(chain_include) > 0:
            list_chains = chain_include.upper().split(',')
        else:
            list_chains = self.config_chains.upper().split(',')

        if len(distance) > 0:
            shell_dist = float(distance)
        else:
            shell_dist = float(self.config_distance)

        if len(key_pdb) > 0:
            pdb_filename = key_pdb
        else:
            pdb_filename = self.pdb_file


        #Create a list of xyz coords that we want to match to (i.e. flatten the pdb dict)
        list_coords = []

        #Loop our chains
        for chain in list_chains:

            #The dict pdb[atom][chain][resi_num] is a list of dicts
            #Only select the atoms and chains we actually want
            for resi_num in dict_pdb[pdb_filename]['dict_atom'][chain]:

                #Loop the resi_nums
                for atom in dict_pdb[pdb_filename]['dict_atom'][chain][resi_num]:

                    #Skip atoms we do not want
                    if atom['atom_name'] not in list_atoms:
                        continue

                    #Else append the x,y,z as a list of lists
                    list_coords.append([atom['x_coor'], atom['y_coor'], atom['z_coor']])

        #Create a dict to work into
        dict_contact = {}

        #Loop our chains
        for chain in list_chains:

            #Check if chain exists in our output dict
            if chain not in dict_contact:
                dict_contact[chain] = {}

            #Loop our dict_pdb
            for resi_num in dict_pdb[pdb_filename]['dict_atom'][chain]:

                #Add the residue to the dict
                if resi_num not in dict_contact[chain]:
                    dict_contact[chain][resi_num] = []

                #Loop each search residue
                for atom in dict_pdb[pdb_filename]['dict_atom'][chain][resi_num]:

                    #Exclude if not in our list
                    if atom['atom_name'] not in list_atoms:
                        continue

                    #Search each residue
                    for search_atm in list_coords:
                        
                        #Exclude if the same
                        if (search_atm[0] == atom['x_coor'] and
                            search_atm[1] == atom['y_coor'] and
                            search_atm[2] == atom['z_coor']):
                            continue

                        #Calculate the Euclidean distance
                        dist = euc_dist(atom['x_coor'], search_atm[0],
                                        atom['y_coor'], search_atm[1],
                                        atom['z_coor'], search_atm[2])

                        #Append to our return dict if dist is below our threshold
                        if dist <= shell_dist:
                            dict_contact[chain][resi_num].append(dist)

        return dict_contact

    def classified_count(self, dict_contact, dict_classified):
        """Count our classifiers"""

        #Import our counter
        from collections import Counter

        """
        Distance
        <= 16
        17 to 24A
        >= 25
        """

        #Setup our lists
        list_16fewer = []
        list_17to24 = []
        list_25greater = []

        #Loop the locations
        for loc in dict_classified:

            #Loop the mutations
            for mut in dict_classified[loc]:

                #Skip stops
                if mut == "*":
                    continue

                #Skip residues without location data
                if loc not in dict_contact[self.report_chain]:
                    continue

                #<= 16
                if len(dict_contact[self.report_chain][loc]) <= 16:
                    list_16fewer.append(dict_classified[self.report_chain][loc][mut])

                #17 to 24
                if len(dict_contact[self.report_chain][loc]) >= 17 and len(dict_contact[self.report_chain][loc]) <= 24:
                    list_17to24.append(dict_classified[loc][mut])

                #>= 25
                if len(dict_contact[self.report_chain][loc]) >= 25:
                    list_25greater.append(dict_classified[loc][mut])

        #Count the lists
        str_return = '\n'.join(map(str, [
            "Contact Number: <= 16",
            pretty_counter_dicts(dict(Counter(list_16fewer))),
            "",
            "Contact Number: 17 to 24",
            pretty_counter_dicts(dict(Counter(list_17to24))),
            "",
            "Contact Number: >= 25",
            pretty_counter_dicts(dict(Counter(list_25greater)))
            ]))

        print(str_return)
        return str_return
   
if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[Contact Number Error] This classifier needs to be ran within the context of PACT.")   
