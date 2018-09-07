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

"""interface distance - calculate the distance to the nearest neighbor at the interface"""

#Set the name of this module
module_name = "Interface Distance"

#Check to see if our version is supported
from sys import version_info
if version_info < (3,4):
    print("[" + module_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError

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

class interface_distance:
    """Calculate the distince to a interface"""

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
            print("[Interface Distance] Using config file options")
            self.config_file_check = False

            try:
                #Set the pdb file name to look for in the pdb file dict
                self.pdb_file = self.config_file.get('interface_distance', 'pdb_file')

                #Select the main chain
                self.main_chain = self.config_file.get('interface_distance', 'main_chain').upper()

                #Select the secondary chain
                self.secondary_chains = [chain.strip(' ') for chain in
                               self.config_file.get('interface_distance', 'secondary_chains').upper().split(',')]

            except NoOptionError:
                print("[Interface Distance Warning] The config file is missing some options")
        else:
            print("[Interface Distance] Not using config file options")
            self.config_file_check = True

        return

    def manual_config(self, dict_configopts):
        """Set manual config options"""

        #Set the pdb file name to look for in the pdb file dict
        self.pdb_file = dict_configopts['pdb_file']

        #Select the main chain
        self.main_chain = dict_configopts['main_chain'].upper()

        #Select the main chain
        self.secondary_chains = [chain.strip(' ') for chain in
                               dict_configopts['secondary_chain'].upper().split(',')]
        
        return

    def secondary_chain(self, dict_pdb):
        """Return a list of contact coords"""

        #Create a list of xyz coords that we want to match to (i.e. flatten the pdb dict)
        return list([atom['x_coor'], atom['y_coor'], atom['z_coor']]
                           for chain in self.secondary_chains
                           for resi_num in dict_pdb[self.pdb_file]['dict_atom'][chain]
                           for atom in dict_pdb[self.pdb_file]['dict_atom'][chain][resi_num])

    def interface_distance(self, dict_pdb, dict_configopts = None):
        """Calculate the contact number of a residue"""

        #Set manual config options
        if dict_configopts != None:
            self.manual_config(dict_configopts)
        elif self.config_file_check and dict_configopts == None:
            print("[" + module_name + " Error] The config file or manual config is missing")
            quit()

        #Return the list of potential contact coords
        nparr_secondary_chain = np.array(self.secondary_chain(dict_pdb))

        #Create a dict to work into
        dict_interface = {self.main_chain:{}}

        #Loop the residues
        for resi_num in dict_pdb[self.pdb_file]['dict_atom'][self.main_chain]:

            #Create two lists of mainchain and sidechain atoms
            list_mainchain = []
            list_sidechain = []

            #Return the mainchain atoms
            for atom in dict_pdb[self.pdb_file]['dict_atom'][self.main_chain][resi_num]:

                #Check for mainchain
                if atom['atom_name'] in ['N', 'CA', 'C', 'O', 'N', 'H']:
                    list_mainchain.append(np.min(cdist(np.array([
                        [atom['x_coor'], atom['y_coor'], atom['z_coor']]]), 
                        nparr_secondary_chain, metric='euclidean')))
                else:
                    list_sidechain.append(np.min(cdist(np.array([
                        [atom['x_coor'], atom['y_coor'], atom['z_coor']]]), 
                        nparr_secondary_chain, metric='euclidean')))

            #Convert to array, negative values go to zero
            nparr_mainchain = np.array(list_mainchain).clip(min=0)
            nparr_sidechain = np.array(list_sidechain).clip(min=0)

            #Calculate the euc dist and report
            dict_interface[self.main_chain][resi_num] = {
                'main_min':np.min(nparr_mainchain),
                'main_mean':np.mean(nparr_mainchain),
                'main_max':np.max(nparr_mainchain),
                'side_min':np.min(nparr_sidechain),
                'side_mean':np.mean(nparr_sidechain),
                'side_max':np.max(nparr_sidechain),
                }

        return dict_interface
  
if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[" + module_name + " Error] This classifier needs to be ran within the context of PACT.")   
