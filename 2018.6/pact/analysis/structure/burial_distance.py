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

"""burial distance - calculate the burial distance of any atom"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Burial Distance Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from math import sqrt, pow, pi, cos, sin
from random import random
from pact.pact_common import pretty_counter_dicts, euc_dist

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class burial_distance:
    """Calculate the exposed residues and then the minimum burial distance to any atom"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Setup the globals"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Get any config file entries
        if self.config_file.has_section("burial_distance"):
            try:
                #Classifier specific variables
                self.pdb_file = self.config_file.get('burial_distance', 'pdb_file')
                self.config_chains = self.config_file.get('burial_distance', 'chains')
                self.report_chain = self.config_file.get('burial_distance', 'report_chain')
                self.num_points = int(self.config_file.get('burial_distance', 'points'))
            except NoOptionError:
                print("[Burial Distance Error] the config file is missing any of the options: pdb_file, atoms, or chains.")
        else:
            print("[Burial Distance Error] the config file is missing the section [burial_distance]")

        return

    def vdw_water_radius(self, atom):
        """Calculate the vdw+h2o radius"""

        #Check for a letter or number at the first position of the atom
        if atom[0].isalpha():
            atom = atom[0]
        else:
            atom = atom[1]

        #Get the van der Walls radius, and water radius
        dict_rvdw = {
            'N':1.4,
            'C':1.548,
            'S':1.808,
            'H':1.1,
            'O':1.348,
            'P':1.88,
            }

        #Calculate the radius (rvdw + rwater)
        if len(atom) > 0:
            #Use the atom from the dict (use the first letter of the atom field as it could be ND1 or other combination)
            return dict_rvdw[atom] + 1.4
        else:
            #Default to the carbon radius given it's more likely
            return 1.548 + 1.4

    def fibonacci_sphere(self, xatm, yatm, zatm, atom = ''):
        """Use the Fibonacci sphere algorithm to make random points around a coordinate"""

        #Get the radius
        radius = self.vdw_water_radius(atom)

        #Setup the variables
        rnd = random() * self.num_points

        list_points = []
        offset = 2./self.num_points
        increment = pi * (3. - sqrt(5.));

        #Loop and create the point
        for i in range(self.num_points):
            y = ((i * offset) - 1) + (offset / 2);
            r = sqrt(1 - pow(y,2))

            phi = ((i + rnd) % self.num_points) * increment

            x = cos(phi) * r
            z = sin(phi) * r

            #Multiple the xyz by the wanted radius then add it to the atom coordinate
            list_points.append([          
                (x * radius) + xatm,
                (y * radius) + yatm,
                (z * radius) + zatm,
                atom
                ])

        return list_points

    def sphere_points(self, dict_pdb, file_pdb, list_chains):
        """Calculate all of the points around the sidechains"""

        #Compute a dict that has [[xyz], chain, loc]
        #of 400 points of all sidechain residues CA to whatever
        #Loop the chains that we want
        print("[Burial Distance] Calculating the points around each atom.")
        dict_points = {}
        for chain in list_chains:

            #Add the chain
            if chain not in dict_points:
                dict_points[chain] = {}

            #Loop the locations
            for resi_num in dict_pdb[file_pdb]['dict_atom'][chain]:

                #Add the residue
                if resi_num not in dict_points[chain]:
                    dict_points[chain][resi_num] = []

                #Loop the resi_nums
                for atom in dict_pdb[file_pdb]['dict_atom'][chain][resi_num]:

                    #Assign our atom name to a temp var
                    atom_name = atom['atom_name']

                    #Check if our atom has a number in the first position use the rest
                    if not atom_name[0].isalpha():
                        atom_name = atom_name[1:]

                    #Append our list
                    dict_points[chain][resi_num] += self.fibonacci_sphere(
                        atom['x_coor'],
                        atom['y_coor'],
                        atom['z_coor'],
                        atom_name                       
                        )

        return dict_points

    def nearest_neighbors(self, dict_pdb, file_pdb, list_chains):
        """Return a dict of [chain][loc][[chain, loc],[]] of residues within 10A"""
        
        #Create a dict to work out of
        list_coords = []

        #Loop the chains
        for chain in list_chains:

            #Loop the locations
            for resi_num in dict_pdb[file_pdb]['dict_atom'][chain]:

                #Loop the resi_nums
                for atom in dict_pdb[file_pdb]['dict_atom'][chain][resi_num]:

                    #Only look at Cb and HA1 for GLY
                    if atom['atom_name'] != "CB" and atom['res_name'] != "GLY":
                        continue

                    if atom['atom_name'] != "HA" and atom['res_name'] == "GLY":
                        continue

                    #Add the coords to the list
                    list_coords.append([atom['x_coor'], atom['y_coor'], atom['z_coor'], chain, resi_num])

        dict_nn = {}
        #Loop the chains
        for chain in list_chains:

            #Add the chain
            if chain not in dict_nn:
                dict_nn[chain] = {}

            #Loop the locations
            for resi_num in dict_pdb[file_pdb]['dict_atom'][chain]:

                #Add the residue
                if resi_num not in dict_nn[chain]:
                    dict_nn[chain][resi_num] = []

                #Loop the resi_nums
                for atom in dict_pdb[file_pdb]['dict_atom'][chain][resi_num]:

                    #Only look at Cb and HA1 for GLY
                    if atom['atom_name'] != "CB" and atom['res_name'] != "GLY":
                        continue

                    if atom['atom_name'] != "HA1" and atom['res_name'] == "GLY":
                        continue

                    #Loop and add residues under 10A
                    for coords in list_coords:

                        #Check if itself
                        if coords[4] == resi_num:
                            continue

                        #Check dist
                        dist = euc_dist(atom['x_coor'], coords[0],
                                        atom['y_coor'], coords[1],
                                        atom['z_coor'], coords[2])

                        #Evaluate
                        if dist < 10:
                            dict_nn[chain][resi_num].append([coords[3], coords[4]])
       
        return dict_nn

    def remove_nn(self, dict_pdb, file_pdb, dict_nn, dict_points):
        """Remove overlapping near neighbors points"""

        #Loop the chains
        for chain in dict_points:

            #Loop the residues locations
            for resi_loc in dict_points[chain]:

                print("Residue: " + str(resi_loc))

                #Get the list of neighbors ["chain", resi_loc], [...
                list_nn = dict_nn[chain][resi_loc]

                #Loop the NN residues
                for nn in range(0, len(list_nn)):

                    #NN Chain
                    nn_chain = list_nn[nn][0]

                    #NN Location
                    nn_location = list_nn[nn][1]

                    #Get the neighbor points [[x,y,z,atom], [x..
                    len_original_points = len(dict_points[nn_chain][nn_location])

                    #Loop the original coors of the structure
                    for atom in dict_pdb[file_pdb]['dict_atom'][chain][resi_loc]:

                        #Create a list to append to
                        list_points_to_save = []

                        #Calculate the rvdw+rw dist
                        radi_dist = self.vdw_water_radius(atom['atom_name'])

                        #Evaluate the points
                        for point in dict_points[nn_chain][nn_location]:

                            #Calculate the distance from atom to nn point
                            dist = euc_dist(
                                point[0], atom['x_coor'],
                                point[1], atom['y_coor'],
                                point[2], atom['z_coor']
                                )

                            #Evaluate if we're overlapping, save if the distance is over
                            if dist > radi_dist:
                               list_points_to_save.append(point)
                               
                        #Overwrite the nn list
                        dict_points[nn_chain][nn_location] = list_points_to_save

        return dict_points

    def remove_all_sites(self, dict_pdb, file_pdb, dict_points):
        """Do a final pass to remove overlapping points"""

        #Loop the chains
        for chain in dict_points:

            #Loop the residues locations
            for resi_loc in dict_points[chain]:

                print("Residue: " + str(resi_loc))

                #Skip points with zero
                if len(dict_points[chain][resi_loc]) == 0:
                    print("Skipped as it is zero.")
                    continue

                #Loop the original coors of the structure
                for atom in dict_pdb[file_pdb]['dict_atom'][chain][resi_loc]:

                    #Calculate the rvdw+rw dist
                    radi_dist = self.vdw_water_radius(atom['atom_name'])

                    #Evaluate search chains
                    for search_chain in dict_points:

                        #Evaluate all locations
                        for search_loc in dict_points[search_chain]:

                            #Create a list to append to
                            list_points_to_save = []

                            #Evaluate the points
                            for point in dict_points[search_chain][search_loc]:

                                #Calculate the distance from atom to nn point
                                dist = euc_dist(
                                    point[0], atom['x_coor'],
                                    point[1], atom['y_coor'],
                                    point[2], atom['z_coor']
                                    )


                                #Evaluate if we're overlapping, save if the distance is over
                                if dist > radi_dist:
                                    list_points_to_save.append(point)
                               
                            #Overwrite the nn list
                            dict_points[search_chain][search_loc] = list_points_to_save

        return dict_points

    def burial_distance(self, dict_pdb, key_pdb='', chain_include=''):
        """Calculate the contact number of a residue"""

        #Use the config file or the function call
        if len(chain_include) > 0:
            list_chains = chain_include.upper().split(',')
        else:
            list_chains = self.config_chains.upper()

        if len(key_pdb) > 0:
            pdb_filename = key_pdb
        else:
            pdb_filename = self.pdb_file


        from scipy.spatial.distance import cdist
        from numpy import array, matrix


        x = array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4], [5,5,5]])
        y = array([[10, 10, 10], [20, 20, 20], [30, 30, 30], [40, 40, 40]])

        dist = cdist(y, x, metric='euclidean')

        print(dist)

        #Subtract
        z = matrix([[1], 
                   [1], 
                   [1],
                   [1]])

        print(z)

        print(dist - z)


        quit()

        #Calculate the nearest neighbors
        #Return a dict of [chain][loc] = [[chain, loc],[]] of residues within 10A
        print("[Burial Distance] Calculating the nearest neighbors to CB or HB1 in 10A shell")
        dict_nn = self.nearest_neighbors(dict_pdb, pdb_filename, list_chains)

        #Calculate the sphere around sidechain carbons or hydrogens for gly
        #[chain][resi_num] = [x,y,z,atom], [...
        print("[Burial Distance] Fibonacci sphere points around all atoms")
        dict_points = self.sphere_points(dict_pdb, pdb_filename, list_chains)

        #Evaluate the neighbors and the points as a first pass
        #[chain][resi_num] = [x,y,z,atom], [...
        print("[Burial Distance] First pass: removing points that overlap nearest neighbors")
        dict_nonn_points = self.remove_nn(dict_pdb, pdb_filename, dict_nn, dict_points)

        #Perform a second pass and evaluate against all sites
        print("[Burial Distance] Second pass: removing points that overlap all residues")
        dict_all_points = self.remove_all_sites(dict_pdb, pdb_filename, dict_nonn_points)

        #Output our surface residues
        print("Surface Resis")
        for chain in dict_all_points:
            for resi_loc in dict_all_points[chain]:
                print(str(resi_loc) + " " + str(len(dict_all_points[chain][resi_loc])))

        return dict_all_points

if __name__ == '__main__':

    #Remind the user that the classifier needs to be ran within the context of PACT
    print("Burial Distance Error] This classifier needs to be ran within the context of PACT.")   
