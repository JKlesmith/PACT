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

#Set the name of this module
module_name = "Contact Number"

#Check to see if our version is supported
from sys import version_info
if version_info < (3,4):
    print("[" + module_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from collections import Counter
from configparser import NoOptionError
from pact.pact_common import pretty_counter_dicts

try:
    import numpy as np
except ImportError:
    print("[Error] Numpy is not installed and is required.")
    quit()

try:
    from scipy.spatial.distance import cdist
except ImportError:
    print("[Error] SciPy is not installed and is required.")
    quit()

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.12"
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

            #Using config file
            print("[Contact Number] Using config file options")
            self.config_file_check = False

            try:
                #Set the pdb file name to look for in the pdb file dict
                self.pdb_file = self.config_file.get('contact_number', 'pdb_file')

                #Select which atoms to base our calculation off of
                self.atom = self.config_file.get('contact_number', 'atom').upper()

                #Select which chains to calculate for
                self.chains = [chain.strip(' ') for chain in
                               self.config_file.get('contact_number', 'chains').upper().split(',')]

                #Define the distance shell for neighbors
                self.distance = float(self.config_file.get('contact_number', 'distance'))

                #Select the chain to count our classifiers on
                self.classifer_chain = self.config_file.get('contact_number', 'classifer_chain').upper()

            except NoOptionError:
                print("[Contact Number Warning] The config file is missing some options")
        else:
            print("[Contact Number] Not using config file options")
            self.config_file_check = True

        return

    def manual_config(self, dict_configopts):
        """Set manual config options"""

        #Set the pdb file name to look for in the pdb file dict
        self.pdb_file = dict_configopts['pdb_file']

        #Select the atom to base our calculation off of
        self.atom = dict_configopts['atom'].upper()

        #Select which chains to calculate for
        self.chains = [chain.strip(' ') for chain in dict_configopts['chains'].upper().split(',')]

        #Define the distance shell for neighbors
        self.distance = float(dict_configopts['distance'])

        #Select the chain to count our classifiers on
        self.classifer_chain = dict_configopts['classifier_chain'].upper()

        return

    def contact_coords(self, dict_pdb):
        """Return a list of contact coords"""                 

        #Create a list of xyz coords that we want to match to (i.e. flatten the pdb dict)
        return list([atom['x_coor'], atom['y_coor'], atom['z_coor']]
                           for chain in self.chains
                           for resi_num in dict_pdb[self.pdb_file]['dict_atom'][chain]
                           for atom in dict_pdb[self.pdb_file]['dict_atom'][chain][resi_num]
                           if atom['atom_name'] == self.atom)

    def contact_number(self, dict_pdb, dict_configopts = None):
        """Calculate the contact number of a residue"""

        #Set manual config options
        if dict_configopts != None:
            self.manual_config(dict_configopts)
        elif self.config_file_check and dict_configopts == None:
            print("[" + module_name + " Error] The config file or manual config is missing")
            quit()

        #Return the list of potential contact coords
        nparr_contact_coords = np.array(self.contact_coords(dict_pdb))

        #Create a dict to work into
        dict_contact = {}

        #Loop our chains
        for chain in self.chains:

            #Check if chain exists in our output dict
            if chain not in dict_contact:
                dict_contact[chain] = {}

            #Loop our dict_pdb
            for resi_num in dict_pdb[self.pdb_file]['dict_atom'][chain]:

                #Loop each search residue
                nparr_atm_coord = np.array([[atom['x_coor'], atom['y_coor'], atom['z_coor']]
                                  for atom in dict_pdb[self.pdb_file]['dict_atom'][chain][resi_num]
                                  if atom['atom_name'] == self.atom])

                #Calculate the euc dist
                resi_dists = cdist(nparr_atm_coord, nparr_contact_coords, metric='euclidean')

                #Find the number of elements that are less than our shell dist, then subtract 1 for itself
                dict_contact[chain][resi_num] = np.count_nonzero(resi_dists <= self.distance) - 1                                    

        return dict_contact

    def classified_count(self, dict_contact, dict_classified, dict_configopts = None):
        """Count our classifiers"""

        #Set manual config options
        if dict_configopts != None:
            self.manual_config(dict_configopts)
        elif self.config_file_check and dict_configopts == None:
            print("[" + module_name + " Error] The config file or manual config is missing")
            quit()

        #Import our counter
        
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
                if dict_contact[self.report_chain][loc] <= 16:
                    list_16fewer.append(dict_classified[self.report_chain][loc][mut])

                #17 to 24
                if dict_contact[self.report_chain][loc] >= 17 and dict_contact[self.report_chain][loc] <= 24:
                    list_17to24.append(dict_classified[loc][mut])

                #>= 25
                if dict_contact[self.report_chain][loc] >= 25:
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
    print("[" + module_name + " Error] This classifier needs to be ran within the context of PACT.")   
