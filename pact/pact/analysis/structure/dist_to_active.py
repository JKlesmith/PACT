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

"""dist_to_active - calculate the distance from any residue to a ligand in the structure"""

#Set the name of this module
module_name = "Distance to Active Site"

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

class dist_to_active:
    """Class with methods to calculate the distance to ligands"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Setup the globals"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Get any config file entries
        if self.config_file.has_section("distance_to_active"):

            #Using config file
            print("[Distance to Active Site] Using config file options")
            self.config_file_check = False

            #Classifier specific variables
            try:
                #Set the pdb file name to look for in the pdb file dict
                self.pdb_file = self.config_file.get('distance_to_active', 'pdb_file')

                #Select which atoms to base our calculation off of
                self.atoms = [atm.strip(' ') for atm in
                             self.config_file.get('distance_to_active', 'atoms').upper().split(',')]

                #Select which chains to calculate for
                self.chains = [chain.strip(' ') for chain in
                               self.config_file.get('distance_to_active', 'chains').upper().split(',')]

                #Define which ligand type to calculate the distance to
                self.active_type = self.config_file.get('distance_to_active', 'active_type').upper()

                #Define which chain the active site residues or ligands are on
                self.active_chains = self.config_file.get('distance_to_active', 
                                                                'active_chains').upper().split(',')

                #Define which ligands to calculate the distance to
                self.active_ligands = [ligand.strip(' ') for ligand in
                                       self.config_file.get('distance_to_active', 'active_ligands').upper().split(',')]

                #Define which active site residues to calculate the distance to
                self.active_resis = list(map(int, self.config_file.get('distance_to_active', 
                                                                             'active_residues').upper().split(',')))

                #Select the chain to count our classifiers on
                self.classifer_chain = self.config_file.get('distance_to_active', 'classifer_chain').upper()

            except NoOptionError:
                print("[Distance to Active Site Warning] The config file is missing some options")

        else:
            #Using config file
            print("[Distance to Active Site] Not using config file options")
            self.config_file_check = True

        return

    def manual_config(self, dict_configopts):
        """Set manual config options"""

        #Set the pdb file name to look for in the pdb file dict
        self.pdb_file = dict_configopts['pdb_file']

        #Select which atoms to base our calculation off of
        self.atoms = [atm.strip(' ') for atm in dict_configopts['atoms'].upper().split(',')]

        #Select which chains to calculate for
        self.chains = [chain.strip(' ') for chain in dict_configopts['chains'].upper().split(',')]

        #Define which ligands to calculate the distance to
        self.active_type = self.config_file.get('distance_to_active', 'active_type').upper()

        #Define which chain the active site residues or ligands are on
        self.active_chains = [chain.strip(' ') for chain in dict_configopts['active_chains'].upper().split(',')]

        #Define which ligands to calculate the distance to
        self.active_ligands = [ligands.strip(' ') for ligands in dict_configopts['active_ligands'].upper().split(',')]

        #Define which active site residues to calculate the distance to
        self.active_resis = list(map(int, dict_configopts['active_residues'].upper().split(',')))

        #Select the chain to count our classifiers on
        self.classifer_chain = dict_configopts['classifier_chain'].upper()

        return

    def active_coords(self, dict_pdb):
        """Return a list of lists of xyz coords for active ligands or residues"""

        #Create a list of lists of hetatm x,y,z values (lists keep order)
        list_active_coord = []

        #The dict pdb[atom][resi_num] is a list of dicts
        if self.active_type == "LIGANDS":

            #Parse the hetatm chains
            for ligand_chain in self.active_chains:

                #Check to see if the ligand chain is in the dict
                if ligand_chain not in dict_pdb[self.pdb_file]['dict_hetatm']:
                    print("[Distance to Active Error] Ligand is not in the pdb.")
                    quit()

                #Parse the hetatms resi_nums
                for hetatm_num in dict_pdb[self.pdb_file]['dict_hetatm'][ligand_chain]:

                    #Parse each atom
                    for hetatm in dict_pdb[self.pdb_file]['dict_hetatm'][ligand_chain][hetatm_num]:

                        #Skip ligands we do not want using the first entry
                        if (hetatm['res_name'] not in self.active_ligands or
                            hetatm['chain'] not in self.active_chains):
                            continue

                        #Add the coords to the list
                        list_active_coord.append([
                            hetatm['x_coor'],
                            hetatm['y_coor'],
                            hetatm['z_coor']
                            ])
        else:
            #Parse the atom chains
            for ligand_chain in self.active_chains:

                #Parse the hetatms resi_nums
                for atom_num in dict_pdb[self.pdb_file]['dict_atom'][ligand_chain]:

                    #Parse each atom
                    for atm in dict_pdb[self.pdb_file]['dict_atom'][ligand_chain][atom_num]:

                        #Skip ligands we do not want using the first entry
                        if (atm['resi_num'] not in self.active_resis or
                            atm['chain'] not in self.active_chains):
                            continue

                        #Add the coords to the list
                        list_active_coord.append([
                            atm['x_coor'],
                            atm['y_coor'],
                            atm['z_coor']
                            ])

        return list_active_coord

    def dta_dist(self, dict_pdb, dict_configopts = None):
        """Calculates the Euclidean distance from each atom to hetatm"""

        #Set manual config options
        if dict_configopts != None:
            self.manual_config(dict_configopts)

        #Get the active ligands/resis coords
        nparr_active_coord = np.array(self.active_coords(dict_pdb))

        #Create a dict to work into
        dict_dtoa_dist = {}

        #Parse the atm chains
        for atm_chain in self.chains:

            #Add the chain
            if atm_chain not in dict_dtoa_dist:
                dict_dtoa_dist[atm_chain] = {}

            #Parse the residue numbers
            for resi_num in dict_pdb[self.pdb_file]['dict_atom'][atm_chain]:

                #The dict pdb[atom][resi_num] is a list of dicts
                #Make a list of atom coords in our resi number
                #Only keep atoms we want and chains we want
                nparr_atm_coord = np.array([[atom['x_coor'], atom['y_coor'], atom['z_coor']]
                                  for atom in dict_pdb[self.pdb_file]['dict_atom'][atm_chain][resi_num]
                                  if atom['atom_name'] in self.atoms and
                                  atom['chain'] in self.chains])
                
                #Calculate the euc dist and return the min per location
                dict_dtoa_dist[atm_chain][resi_num] = float(np.min(cdist(
                    nparr_atm_coord, nparr_active_coord, metric='euclidean')))

        #Return a dict of {chain:resi_num:min dist]
        return dict_dtoa_dist

    def classified_count(self, dict_dtoa, dict_classified, dict_configopts = None):
        """Count our classifiers"""

        #Set manual config options
        if dict_configopts != None:
            self.manual_config(dict_configopts)
        elif self.config_file_check and dict_configopts == None:
            print("[" + module_name + " Error] The config file or manual config is missing")
            quit()       

        """
        Distance
        < 10A
        >= 10 to 14.999A
        >= 15 to 19.999A
        >= 20A
        """

        #Setup our lists
        list_10closer = []
        list_10to14 = []
        list_15to19 = []
        list_20further = []

        #Loop the locations
        for loc in dict_classified:

            #Loop the mutations
            for mut in dict_classified[loc]:

                #Skip stops
                if mut == "*":
                    continue

                #Skip residues without location data
                if loc not in dict_dtoa[self.classifer_chain]:
                    continue

                #< 10A
                if dict_dtoa[self.classifer_chain][loc] < 10:
                    list_10closer.append(dict_classified[loc][mut])

                #>= 10 to 14.999A
                if dict_dtoa[self.classifer_chain][loc] >= 10 and dict_dtoa[self.classifer_chain][loc] < 15:
                    list_10to14.append(dict_classified[loc][mut])

                #>= 15 to 19.999A
                if dict_dtoa[self.classifer_chain][loc] >= 15 and dict_dtoa[self.classifer_chain][loc] < 20:
                    list_15to19.append(dict_classified[loc][mut])

                #>= 20A
                if dict_dtoa[self.classifer_chain][loc] >= 20:
                    list_20further.append(dict_classified[loc][mut])

        #Count the lists
        str_return = '\n'.join(map(str, [
            "Distance: < 10A",
            pretty_counter_dicts(dict(Counter(list_10closer))),
            "",
            "Distance: >= 10 to 14.999A",
            pretty_counter_dicts(dict(Counter(list_10to14))),
            "",
            "Distance: >= 15 to 19.999A",
            pretty_counter_dicts(dict(Counter(list_15to19))),
            "",
            "Distance: >= 20A",
            pretty_counter_dicts(dict(Counter(list_20further)))
            ]))

        print(str_return)
        return str_return

if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[" + module_name + " Error] This classifier needs to be ran within the context of PACT.")   
