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

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Distance to Active Site Error] Your Python interpreter is too old."
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
            try:
                #Classifier specific variables
                self.pdb_file = self.config_file.get('distance_to_active', 'pdb_file')
                self.config_atoms = self.config_file.get('distance_to_active', 'atoms')
                self.config_chains = self.config_file.get('distance_to_active', 'chains')
                self.report_chain = self.config_file.get('distance_to_active', 'report_chain')
                self.config_ligands = self.config_file.get('distance_to_active', 'ligands')
                self.config_activeresis = self.config_file.get('distance_to_active', 'active_residues')
                self.config_ligandchains = self.config_file.get('distance_to_active', 'ligand_chains')
            except NoOptionError:
                print("[Distance to Active Site Error] the config file is missing any of the options: pdb_file, atoms, chains, or ligands.")
        else:
            print("[Distance to Active Site Error] the config file is missing the section [distance_to_active]")

        return

    def dta_dist(self, dict_pdb, key_pdb='', atom_include='', chain_include='', ligand_include='', ligand_chains='', active_residues=''):
        """Calculates the Euclidean distance from each atom to hetatm"""

        #Use the config file or the function call
        if len(atom_include) > 0:
            list_atoms = atom_include.upper().split(',')
        else:
            list_atoms = self.config_atoms.upper().split(',')

        if len(chain_include) > 0:
            list_chains = chain_include.upper().split(',')
        else:
            list_chains = self.config_chains.upper().split(',')

        if len(ligand_include) > 0:
            list_ligands = ligand_include.upper().split(',')
        else:
            list_ligands = self.config_ligands.upper().split(',')

        if len(active_residues) > 0:
            list_residues = list(map(int, active_residues.upper().split(',')))
        else:
            if len(self.config_activeresis) > 0:
                list_residues = list(map(int, self.config_activeresis.split(',')))
            else:
                list_residues = None

        if len(key_pdb) > 0:
            pdb_filename = key_pdb
        else:
            pdb_filename = self.pdb_file

        if len(ligand_chains) > 0:
            list_ligandchains = ligand_chains.upper().split(',')
        else:
            list_ligandchains = self.config_ligandchains.upper().split(',')

        #Create a list of lists of hetatm x,y,z values (lists keep order)
        list_active_coord = []

        #The dict pdb[atom][resi_num] is a list of dicts
        if len(list_ligands[0]) > 0:
            #Parse the hetatm chains
            for ligand_chain in list_ligandchains:

                #Parse the hetatms resi_nums
                for hetatm_num in dict_pdb[pdb_filename]['dict_hetatm'][ligand_chain]:

                    #Parse each atom
                    for hetatm in dict_pdb[pdb_filename]['dict_hetatm'][ligand_chain][hetatm_num]:

                        #Skip ligands we do not want using the first entry
                        if (hetatm['res_name'] not in list_ligands or
                            hetatm['chain'] not in list_ligandchains):
                            continue

                        #Add the coords to the list
                        list_active_coord.append([
                            hetatm['x_coor'],
                            hetatm['y_coor'],
                            hetatm['z_coor']
                            ])
        else:
            #Parse the atom chains
            for ligand_chain in list_ligandchains:

                #Parse the hetatms resi_nums
                for atom_num in dict_pdb[pdb_filename]['dict_atom'][ligand_chain]:

                    #Parse each atom
                    for atm in dict_pdb[pdb_filename]['dict_atom'][ligand_chain][atom_num]:

                        #Skip ligands we do not want using the first entry
                        if (atm['resi_num'] not in list_residues or
                            atm['chain'] not in list_ligandchains):
                            continue

                        #Add the coords to the list
                        list_active_coord.append([
                            atm['x_coor'],
                            atm['y_coor'],
                            atm['z_coor']
                            ])

        #Create a dict to work into
        dict_dtoa_dist = {}

        #Parse the atm chains
        for atm_chain in list_chains:

            #Add the chain
            if atm_chain not in dict_dtoa_dist:
                dict_dtoa_dist[atm_chain] = {}

            #Parse the residue numbers
            for resi_num in dict_pdb[pdb_filename]['dict_atom'][atm_chain]:

                #Add the residue to the dict
                if resi_num not in dict_dtoa_dist[atm_chain]:
                    dict_dtoa_dist[atm_chain][resi_num] = []

                #The dict pdb[atom][resi_num] is a list of dicts
                #Let's make a new list with the Euclidean distance of the atoms from list atm to atms in list hetatm
                for atom in dict_pdb[pdb_filename]['dict_atom'][atm_chain][resi_num]:

                    #Skip atoms we do not want
                    if atom['atom_name'] not in list_atoms:
                        continue

                    #Skip chains we do not want
                    if atom['chain'] not in list_chains:
                        continue

                    #Loop each hetatm
                    for hetatm in list_active_coord:
               
                        #Calculate the Euclidean distance
                        dist = euc_dist(atom['x_coor'], hetatm[0],
                                        atom['y_coor'], hetatm[1],
                                        atom['z_coor'], hetatm[2])

                        dict_dtoa_dist[atm_chain][resi_num].append(dist)

        #Return a dict of {resi_num:[dist, dist, dist]...
        return dict_dtoa_dist

    def classified_count(self, dict_dtoa, dict_classified):
        """Count our classifiers"""

        #Import our counter
        from collections import Counter

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
                if loc not in dict_dtoa[self.report_chain]:
                    continue

                #< 10A
                if min(dict_dtoa[self.report_chain][loc]) < 10:
                    list_10closer.append(dict_classified[loc][mut])

                #>= 10 to 14.999A
                if min(dict_dtoa[self.report_chain][loc]) >= 10 and min(dict_dtoa[self.report_chain][loc]) < 15:
                    list_10to14.append(dict_classified[loc][mut])

                #>= 15 to 19.999A
                if min(dict_dtoa[self.report_chain][loc]) >= 15 and min(dict_dtoa[self.report_chain][loc]) < 20:
                    list_15to19.append(dict_classified[loc][mut])

                #>= 20A
                if min(dict_dtoa[self.report_chain][loc]) >= 20:
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
    print("[Distance to Active Site Error] This classifier needs to be ran within the context of PACT.")   
