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

#FUTURE: DSSP will fail at BASP, BSER (4 char residue codes), pre-filter

"""pdb_import - import PDB files for use in other scripts"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[PDB Import Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError, NoSectionError
from pact.pact_common import file_checker

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class pdb_import:

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("pdb_import"):
            print("[PDB Import Error] the config file is missing the section for PDB Import")

        #Set the number of pdb files
        try:
            self.numfiles = int(self.config_file.get('pdb_import', 'numpdb'))
        except ValueError:
            print("[PDB Import Error] numpdb is not set properly in the config file.")
        except NoOptionError:
            print("[PDB Import Error] The pdb_import config section is incorrect.")
            print("[PDB Import Error] Missing: numpdb")
            quit()

        #Get the pdb names
        self.dict_pdbs = {}

        #Validate each pdb file and store the full filename
        try:
            #Loop each file
            for i in range(1, self.numfiles + 1):
                #Get the pdb filename
                file_pdb = self.config_file.get('pdb_import', 'file_' + str(i))

                #See if the file exists
                if file_checker(self.dict_protocolconfig['directory'] + file_pdb):
                    #Create a dict to work into
                    self.dict_pdbs[file_pdb] = {}

                    #Store the valid filename
                    self.dict_pdbs[file_pdb] = {'filename':self.dict_protocolconfig['directory'] + file_pdb}
        except NoOptionError:
            print("[PDB Import Error] The pdb_import config file is incorrect.")
            print("[PDB Import Error] Missing: file_#")
            quit()
        
        return

    def pdb_atoms(self, file_pdb):
        """Load the PDB atoms into two dicts"""

        #Open the file
        with open(self.dict_pdbs[file_pdb]['filename'], 'r') as infile:
            list_lines = infile.read().splitlines()
        
        #Create the dicts for the atom and hetatms
        dict_atom = {}
        dict_hetatom = {}

        #Iterate through each line
        for line in list_lines:

            #Handle Atom vs Hetatm
            if line.startswith("ATOM"):

                #Get the chain
                chain = line[21:22].strip()

                #Check if the chain is in the dict
                if chain not in dict_atom:
                    dict_atom[chain] = {}

                #Get the residue number
                resi_num = int(line[22:26].strip())

                #Add the residue to the dict
                if resi_num not in dict_atom[chain]:
                    dict_atom[chain][resi_num] = []

                #Add the information to the residue
                dict_atom[chain][resi_num].append({
                'atom_id':int(line[6:11].strip()),
                'atom_name':line[12:16].strip(),
                'alt_loc':line[16:17].strip(),
                'res_name':line[17:20].strip(),
                'chain':line[21:22].strip(),
                'resi_num':resi_num,
                'x_coor':float(line[30:38].strip()),
                'y_coor':float(line[38:46].strip()),
                'z_coor':float(line[46:54].strip()),
                'occupancy':float(line[54:60].strip()),
                'temp':float(line[60:66].strip()),
                'element':line[76:78].strip(),
                'charge':line[78:80].strip()})

            elif line.startswith("HETATM"):

                #Get the chain
                chain = line[21:22].strip()

                #Check if the chain is in the dict
                if chain not in dict_hetatom:
                    dict_hetatom[chain] = {}

                #Get the residue number
                resi_num = int(line[22:26].strip())

                #Add the residue to the dict
                if resi_num not in dict_hetatom[chain]:
                    dict_hetatom[chain][resi_num] = []

                #Add the information to the residue
                dict_hetatom[chain][resi_num].append({
                'atom_id':int(line[6:11].strip()),
                'atom_name':line[12:16].strip(),
                'alt_loc':line[16:17].strip(),
                'res_name':line[17:20].strip(),
                'chain':line[21:22].strip(),
                'resi_num':resi_num,
                'x_coor':float(line[30:38].strip()),
                'y_coor':float(line[38:46].strip()),
                'z_coor':float(line[46:54].strip()),
                'occupancy':float(line[54:60].strip()),
                'temp':float(line[60:66].strip()),
                'element':line[76:78].strip(),
                'charge':line[78:80].strip()})

        return dict_atom, dict_hetatom

    def dssp(self, file_pdb):
        """Interface with the DSSP program, load all information into a dict"""

        #Import the subprocess module
        from subprocess import check_output

        #Call DSSP, decode the bytes, and split on newline
        dssp_output = check_output([self.dict_programs['dssp'], self.dict_pdbs[file_pdb]['filename']]).decode().split('\r\n')

        #Create a dict to return to
        dict_dssp = {}

        #Find when to start parsing the output
        for i in range(0, len(dssp_output)):
            if dssp_output[i].startswith("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC"):
                startline = i
                break

        #Loop each residue
        for line in dssp_output[startline + 1:-1]:

            #Check for chain breaks (TERS)
            if line[13:14] == "!":
                continue

            #Get the residue number
            try:
                num_residue = int(line[5:10].strip())
            except ValueError:
                print("Error at line:")
                print(line)
                quit()

            #Get the chain letter
            try:
                chain = line[11:12].strip().upper()
            except ValueError:
                print("Error at line:")
                print(line)
                quit()

            #Make a dict for that chain if not present
            if chain not in dict_dssp:
                dict_dssp[chain] = {}

            #Make a dict for that location if not present
            if num_residue not in dict_dssp[chain]:
                dict_dssp[chain][num_residue] = {}

            #Import the dssp information
            dict_dssp[chain][num_residue] = {
            'location':num_residue,
            'chain':line[11:12].strip(),
            'residue':line[13:14].strip(),
            'secondary':line[16:17].strip(),
            'acc':int(line[35:38].strip()),
            'rasa':self.rasa(line[13:14].strip(), int(line[35:38].strip())),
            'frac_burial':self.frac_burial(self.rasa(line[13:14].strip(), int(line[35:38].strip())))
            }

        return dict_dssp

    def rasa(self, residue, acc):
        """For a given acc, return a rasa via the correction by Tien MZ ... Wilke CO (2013)"""
        
        #Setup the Gly-X-Gly table
        rasa_emperical = {
            "A":121, "R":265, "N":187, "D":187, "C":148,
            "E":214, "Q":214, "G":97,  "H":216, "I":195,
            "L":191, "K":230, "M":203, "F":228, "P":154,
            "S":143, "T":163, "W":264, "Y":255, "V":165}

        return acc / rasa_emperical[residue]

    def frac_burial(self, rasa):
        """For a given RASA convert it to fraction buried"""

        #Return the RASA, correct if negative
        if (1 - rasa) < 0:
            return 0
        else:
            return 1 - rasa

    def pdb_import(self):
        """The main function to make a dict of pdb information"""

        #Loop each pdb file
        for pdb_file in self.dict_pdbs:

            #Init our three dicts
            self.dict_pdbs[pdb_file]['dict_atom'] = {}
            self.dict_pdbs[pdb_file]['dict_hetatm'] = {}
            self.dict_pdbs[pdb_file]['dssp'] = {}

            #Get the atom dicts
            dict_atom, dict_hetatm = self.pdb_atoms(pdb_file)

            self.dict_pdbs[pdb_file]['dict_atom'] = dict_atom
            self.dict_pdbs[pdb_file]['dict_hetatm'] = dict_hetatm

            #Get the dssp dicts
            self.dict_pdbs[pdb_file]['dssp'] = self.dssp(pdb_file)

        return self.dict_pdbs

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[PDB Import Error] This script needs to be ran within the context of PACT.")
