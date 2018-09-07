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

#Set the name of this module
module_name = "Burial Distance"

#Check to see if our version is supported
from sys import version_info
if version_info < (3,4):
    print("[" + module_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from math import sqrt, pi
from random import random

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

            #Using config file
            print("[Burial Distance] Using config file options")
            self.config_file_check = False

            try:
                #Set the pdb file name to look for in the pdb file dict
                self.pdb_file = self.config_file.get('burial_distance', 'pdb_file')

                #Select which chains to calculate for
                self.chains = [chain.strip(' ') for chain in
                               self.config_file.get('burial_distance', 'chains').upper().split(',')]

                #Define the number of points to use
                self.num_points = int(self.config_file.get('burial_distance', 'num_points'))

                #Select the chain to count our classifiers on
                self.classifer_chain = self.config_file.get('burial_distance', 'classifer_chain').upper()

            except NoOptionError:
                print("[Burial Distance Warning] The config file is missing some options")
        else:
            print("[Burial Distance] Not using config file options")
            self.config_file_check = True
        
        return

    def manual_config(self, dict_configopts):
        """Set manual config options"""

        #Set the pdb file name to look for in the pdb file dict
        self.pdb_file = dict_configopts['pdb_file']

        #Select which chains to calculate for
        self.chains = [chain.strip(' ') for chain in dict_configopts['chains'].upper().split(',')]

        #Define the number of points to use
        self.num_points = int(dict_configopts['num_points'])

        #Select the chain to count our classifiers on
        self.classifer_chain = dict_configopts['classifier_chain'].upper()

        return

    def atoms_dict_to_list(self, dict_pdb):
        """Return a 2D list of all atoms and atom types"""

        #Create a list of xyz coords that we want to match to (i.e. flatten the pdb dict)
        list_atomtype = []
        list_atoms = []

        #Loop our chains
        for chain in self.chains:

            #The dict pdb[atom][chain][resi_num] is a list of dicts
            #Only select the atoms and chains we actually want
            for resi_num in dict_pdb[self.pdb_file]['dict_atom'][chain]:

                #Loop the atoms
                for atom in dict_pdb[self.pdb_file]['dict_atom'][chain][resi_num]:

                    list_atoms.append([atom['x_coor'], atom['y_coor'], atom['z_coor']])
                    list_atomtype.append(atom['atom_name'])

        return list_atomtype, list_atoms

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

    def fibonacci_sphere(self, list_atoms, list_atomtype):
        """Use the Fibonacci sphere algorithm to make random points around a coordinate"""

        offset = 2/self.num_points
        increment = pi * (3 - sqrt(5))

        #Make a random seed
        rnd = random() * self.num_points

        #Preallocate a numpy array counting 0 to num_points in 1's (array of integers)
        #for i in range(num_points):
        nparr = np.arange(0, self.num_points)

        #   y = ((i * offset) - 1) + (offset / 2)
        nparr_y = np.add(np.subtract(np.multiply(nparr, offset), 1), offset / 2)

        #   r = sqrt(1 - pow(y,2))
        nparr_r = np.sqrt(np.subtract(1, np.power(nparr_y,2)))

        #   phi = ((i + rnd) % num_points) * increment
        nparr_phi = np.multiply(np.mod(np.add(nparr, rnd), self.num_points), increment)

        #   x = cos(phi) * r
        nparr_x = np.multiply(np.cos(nparr_phi), nparr_r)

        #   z = sin(phi) * r
        nparr_z = np.multiply(np.sin(nparr_phi), nparr_r)


        #Append our VDW+Water radius to a list (1D array) for each atom
        list_vdwdist = np.array([self.vdw_water_radius(atom) for atom in list_atomtype])


        #Multiply the XYZ points by the VDW+Water Radius after converting points from 1d to 2d
        #Columns = Different Atoms, Rows = Points, shape = (n atoms, num_points)
        nparr_xrad = np.multiply(nparr_x.reshape((nparr_x.shape[0], 1)), list_vdwdist)
        nparr_yrad = np.multiply(nparr_y.reshape((nparr_y.shape[0], 1)), list_vdwdist)
        nparr_zrad = np.multiply(nparr_z.reshape((nparr_z.shape[0], 1)), list_vdwdist)

        #Add the XYZ points by the original XYZ coords
        #Columns = Different Atoms, Rows = Points, shape = (n atoms, num_points)
        nparr_x_fibvdw = np.add(nparr_xrad, np.array(list_atoms).T[0])
        nparr_y_fibvdw = np.add(nparr_yrad, np.array(list_atoms).T[1])
        nparr_z_fibvdw = np.add(nparr_zrad, np.array(list_atoms).T[2])

        return nparr_x_fibvdw, nparr_y_fibvdw, nparr_z_fibvdw, list_vdwdist

    def surface_points(self, nparr_atomxyz, list_vdwdist, nparr_points):
        """Return an array of surface xyz points"""

        #To save memory we will iterate through each atom position
        #Essentially we will only calculate the distance for non-self positions
        #X-axis is atom positions
        #Y-axis is points
        #Therefore, if our num_points = 100 and position = 1 then ignore column 1 rows 0 to 100
        #Then, subtract the VDW+Water Dist and take np.amin(nparr, axis=1) along x

        #Make a list of the indicies of the aray using the size of nparr xyz
        indices = np.arange(nparr_atomxyz.shape[0])

        #Convert the VDW array from [....] to [[],[]] 1d to 2d (cannot slice as 1d)
        nparr_vdw = list_vdwdist.reshape((list_vdwdist.shape[0], 1))

        #Preallocate a numpy array to copy to
        nparr_bool = np.empty(self.num_points * nparr_atomxyz.shape[0], dtype=bool)

        #Calculate the Euc dist and row min excluding the atom column i and rows point i*num_points
        #This is a slow function, the problem is that 1) cannot broadcast cdist at once
        #or will get a memory error therefore you must window the analysis.
        #The subtraction is also relys upon the index so chaining the index is critical.
        #Expected time = 30sec for 30 points on ~ 6000 atoms
        for i in range(0, nparr_atomxyz.shape[0]):

            try:
                #Calculate the Euclidean distance
                nparr_eucdist = cdist(nparr_points[i * self.num_points : (i * self.num_points) + self.num_points], 
                                      nparr_atomxyz[indices != i, :], metric='euclidean')
            except MemoryError:
                print("[PACT Error] Out of memory, try fewer points.")
                quit()

            #Subtract the VDW+Water distance
            nparr_subdist = np.subtract(nparr_eucdist, nparr_vdw[indices != i, :].T)

            #Return true if point can not be overlapped
            np.copyto(nparr_bool[i * self.num_points : (i * self.num_points) + self.num_points], np.all(nparr_subdist > 0, axis = 1))

        #Return a mask of the original array 
        return nparr_bool

    def burial_distance(self, dict_pdb, dict_configopts = None):
        """Calculate the contact number of a residue"""

        #Set manual config options
        if dict_configopts != None:
            self.manual_config(dict_configopts)
        elif self.config_file_check and dict_configopts == None:
            print("[" + module_name + " Error] The config file or manual config is missing")
            quit()

        #Process the pdb dict to return a list of all atoms and atom types
        print("[Burial Distance] Converting the input PDB")
        list_atomtype, list_atoms = self.atoms_dict_to_list(dict_pdb)

        #Convert the list_atoms to a np array
        nparr_atomxyz = np.array(list_atoms)

        #Make a numpy matrix for points around all atoms
        #Returns xyz arrays of cols = atoms and rows = points
        print("[Burial Distance] Fibonacci sphere points around all atoms")
        nparr_x_fibvdw, nparr_y_fibvdw, nparr_z_fibvdw, list_vdwdist = self.fibonacci_sphere(list_atoms, list_atomtype)
        
        #Convert the matrix into a 2D [x,y,z]
        #Flatten the arrays to 1D using .ravel(), convert into an array, transpose
        #Should be n rows by three columns
        nparr_points = np.array([nparr_x_fibvdw.ravel(), nparr_y_fibvdw.ravel(), nparr_z_fibvdw.ravel()]).T

        #Calculate the Euclidean distance and return a mask for surface points
        print("[Burial Distance] Euclidean distance between points and atoms")
        nparr_surface_mask = self.surface_points(nparr_atomxyz, list_vdwdist, nparr_points)
        nparr_surface_points = nparr_points[nparr_surface_mask]

        print("[Burial Distance] Number of atoms: " + str(nparr_atomxyz.shape[0]))
        print("[Burial Distance] Number of examined points: " + str(nparr_points.shape[0]))
        print("[Burial Distance] Number of surface points: " + str(nparr_surface_points.shape[0]))
        

        #Calculate the fraction of points removed for mainchain and sidechain atoms
        """
        FUTURE: Calculate the fraction of points removed for mainchain and sidechain atoms
        Do what DSSP does.
        points existing / num total points for main or side chain atoms = asa
        Use nparr_surface_mask
        Since the below uses a dict and is not in order must create a new dict matching the input lists
        Could help with setting distance to zero for surface residues instead of 0.0001 cutoff
        """
        
        #Return a dict of [chain][location][sidechain/mainchain _ min/max/mean distance, surface?]
        dict_bd = {}

        #Loop the chains
        for chain in self.chains:

            #Add the chain to our new dict
            if chain not in dict_bd:
                dict_bd[chain] = {}

            #Loop the residues
            for resi_num in dict_pdb[self.pdb_file]['dict_atom'][chain]:

                #Create two lists of mainchain and sidechain atoms
                list_mainchain = []
                list_sidechain = []

                #Return the mainchain atoms
                for atom in dict_pdb[self.pdb_file]['dict_atom'][chain][resi_num]:

                    #Return the vdw+h2o radius
                    vdw_radi = self.vdw_water_radius(atom['atom_name'])

                    #Check for mainchain
                    if atom['atom_name'] in ['N', 'CA', 'C', 'O', 'N', 'H']:
                        list_mainchain.append(np.min(cdist(np.array([
                            [atom['x_coor'], atom['y_coor'], atom['z_coor']]]), 
                            nparr_surface_points, metric='euclidean') - vdw_radi))
                    else:
                        list_sidechain.append(np.min(cdist(np.array([
                            [atom['x_coor'], atom['y_coor'], atom['z_coor']]]), 
                            nparr_surface_points, metric='euclidean') - vdw_radi))

                #Convert to array, negative values go to zero
                nparr_mainchain = np.array(list_mainchain).clip(min=0)
                nparr_sidechain = np.array(list_sidechain).clip(min=0)

                #Calculate the euc dist and report
                #If srfexp dist is < 0.0001 then classify as direct surface exposed
                dict_bd[chain][resi_num] = {
                    'main_min':np.min(nparr_mainchain),
                    'main_mean':np.mean(nparr_mainchain),
                    'main_max':np.max(nparr_mainchain),
                    'side_min':np.min(nparr_sidechain),
                    'side_mean':np.mean(nparr_sidechain),
                    'side_max':np.max(nparr_sidechain),

                    'main_min_surface':np.all(np.min(nparr_mainchain) < 0.0001),
                    'main_mean_surface':np.all(np.mean(nparr_mainchain) < 0.0001),
                    'main_max_surface':np.all(np.max(nparr_mainchain) < 0.0001),
                    'side_min_surface':np.all(np.min(nparr_sidechain) < 0.0001),
                    'side_mean_surface':np.all(np.mean(nparr_sidechain) < 0.0001),
                    'side_max_surface':np.all(np.max(nparr_sidechain) < 0.0001),
                    }

        #Return dict, chain, loc
        return dict_bd

if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[" + module_name + " Error] This classifier needs to be ran within the context of PACT.")   
