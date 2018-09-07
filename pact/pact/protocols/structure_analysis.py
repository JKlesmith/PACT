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

"""structure analysis - basic structure analysis"""

from sys import version_info

#Setup our protocol name
str_protocol_name = "Structure Analysis"

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Protocols:" + str_protocol_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from collections import Counter
from configparser import ConfigParser, NoSectionError, NoOptionError
from os import access, W_OK, R_OK
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class structure_analysis:
    """Basic structure analysis"""

    def __init__(self, file_config, dict_programs, preamble):
        """Assign the class variables from our config file."""

        #Read our config file
        obj_cfgparser = ConfigParser()
        obj_cfgparser.read(file_config)

        #Send the cfgparser object as a global class object
        self.obj_cfgparser = obj_cfgparser

        #Assign our class variables from the config file
        self.output_prefix = obj_cfgparser.get("global", "output_prefix")
        self.directory = obj_cfgparser.get("global", "directory")

        #Read in the workflow
        try:
            self.dict_workflow = {mapping[0].lower(): literal_eval(mapping[1]) 
                              for mapping in obj_cfgparser.items("workflow")}
        except NoSectionError:
            print("[Protocols:" + str_protocol_name + " Error] The config file is incorrect.")
            print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:" + str_protocol_name + " Error] The config file is incorrect.")
            print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [workflow] section (spelling???).")
            quit()

        #Get the PACT preamble
        self.pact_preamble = preamble

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Add trailing slash on directory, windows can be formatted as linux
        if self.directory[-1:] == "\\":
            self.directory = self.directory[:-1] + "/"

        if self.directory[-1:] != "/":
            self.directory = self.directory + "/"

        #Test if directory is writeable
        if not access(self.directory, W_OK) or not access(self.directory, R_OK):
            print("[PACT Error] Output directory is not readable or writable by the program.")
            quit()
        return

    def version(self):
        """Return the version of the protocol"""
        return __version__

    def protocol(self):
        """Main entrypoint for the protocol"""

        #Create a output log file that we can append to
        with open(self.directory + self.output_prefix + "_" + strftime("%m_%d_%Y") + "-" +
                 strftime("%H_%M_%S") + '_output.txt', 'w') as file_output:
            file_output.write(self.pact_preamble + "\n")

            """
            *****************************************
            PDB Import Section (Required)
            *****************************************
            """
            #Check to see if the section is there
            if not self.obj_cfgparser.has_section('pdb_import'):           
                print("[Protocols:" + str_protocol_name + " Error] The pdb_import config file is incorrect.")
                print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [pdb_import] section.")
                quit()

            #Import our combinepact class
            try:
                from pact.analysis.pdb_import import pdb_import
            except ImportError:
                print("[Protocols:" + str_protocol_name + " Error] pdb_import was not found.")

            #Create the object then call the merger
            obj_pdb = pdb_import(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

            #The dict will be like {'pdb name': {data...
            dict_pdb = obj_pdb.pdb_import()

            """
            *****************************************
            Burial Distance (Return the minimum burial distance of that residue)
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['burial_distance']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('burial_distance'):           
                    print("[Protocols:" + str_protocol_name + " Error] The burial_distance config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [burial_distance] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.burial_distance import burial_distance
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.burial_distance was not found.")
                    quit()
        
                #Create our object
                obj_bd = burial_distance(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_bd = obj_bd.burial_distance(dict_pdb)

            """
            *****************************************
            Distance to Active Site
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['distance_to_active']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('distance_to_active'):           
                    print("[Protocols:" + str_protocol_name + " Error] The distance_to_active config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [distance_to_active] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.dist_to_active import dist_to_active
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.dist_to_active was not found.")
                    quit()
        
                #Create our object
                obj_dtoa = dist_to_active(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_dtoa_dist = obj_dtoa.dta_dist(dict_pdb)

            """
            *****************************************
            Contact Number
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['contact_number']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('contact_number'):           
                    print("[Protocols:" + str_protocol_name + " Error] The contact_number config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [contact_number] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.contact_number import contact_number
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.contact_number was not found.")
                    quit()
        
                #Create our object
                obj_contact = contact_number(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_contact = obj_contact.contact_number(dict_pdb) 
                
            """
            *****************************************
            Interface Distance
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['interface_distance']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('interface_distance'):           
                    print("[Protocols:" + str_protocol_name + " Error] The interface_distance config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [interface_distance] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.interface_distance import interface_distance
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.interface_distance was not found.")
                    quit()
        
                #Create our object
                obj_int_dist = interface_distance(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_int_dist = obj_int_dist.interface_distance(dict_pdb)
                
            """
            *****************************************
            Interface Distance
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['contact_map']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('contact_map'):           
                    print("[Protocols:" + str_protocol_name + " Error] The contact_map config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [contact_map] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.contact_map import contact_map
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.contact_map was not found.")
                    quit()
        
                #Create our object
                obj_contact_map = contact_map(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_contact_map = obj_contact_map.contact_map(dict_pdb)

            """
            *****************************************
            HPatch KNN
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['hpatch_knn']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('hpatch_knn'):           
                    print("[Protocols:" + str_protocol_name + " Error] The hpatch_knn config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [hpatch_knn] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.hpatch_knn import hpatch_knn
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.hpatch_knn was not found.")
                    quit()
        
                #Create our object
                obj_hpatchknn = hpatch_knn(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_hpatchknn = obj_hpatchknn.hpatch_knn(dict_pdb) 
            """
            *****************************************
            Dataset Dump CSV
            *****************************************
            """
            #Get the info from the config file
            try:
                chain = self.obj_cfgparser.get('structure_analysis', 'chain').upper()
                pdb_file = self.obj_cfgparser.get('structure_analysis', 'pdb_file')
            except NoOptionError:
                print("[PACT Error] Missing config file option")
                quit()

            #Calculate the number of sites
            list_sites = sorted([x for x in dict_pdb[pdb_file]['dssp'][chain]])
            print("[PDB Info] Found " + str(len(list_sites)) + " residues.")

            #Get the header
            str_output = ','.join([
                "Location",
                "Residue",
                "fraction_burial",
                "dist_to_active",
                "contact_number",
                
                "burial_dist_side_min",
                "burial_dist_side_mean",
                "burial_dist_side_max",
                "burial_dist_main_min",
                "burial_dist_main_mean",
                "burial_dist_main_max",

                "burial_dist_side_min_surf",
                "burial_dist_side_mean_surf",
                "burial_dist_side_max_surf",
                "burial_dist_main_min_surf",
                "burial_dist_main_mean_surf",
                "burial_dist_main_max_surf",

                "interf_dist_side_min",
                "interf_dist_side_mean",
                "interf_dist_side_max",
                "interf_dist_main_min",
                "interf_dist_main_mean",
                "interf_dist_main_max",
                ]) + "\n"

            #Loop the locations
            for loc in list_sites:

                #Get the location
                str_output = str_output + str(loc) + ','

                #Get the residue
                str_output = str_output + dict_pdb[pdb_file]['dssp'][chain][loc]['residue'] + ','

                #Get the fraction burial
                str_output = str_output + str(dict_pdb[pdb_file]['dssp'][chain][loc]['frac_burial']) + ","

                #Dist to active site
                if self.dict_workflow['distance_to_active']:
                    str_output = str_output + str(dict_dtoa_dist[chain][loc]) + ','
                else:
                    str_output = str_output + ','
                    
                #Contact number
                if self.dict_workflow['contact_number']:
                    str_output = str_output + str(dict_contact[chain][loc]) + ','
                else:
                    str_output = str_output + ','

                #Burial Distance
                if self.dict_workflow['burial_distance']:
                    str_output = str_output + ','.join(map(str, [
                        dict_bd[chain][loc]['side_min'],                                   
                        dict_bd[chain][loc]['side_mean'],
                        dict_bd[chain][loc]['side_max'],
                        dict_bd[chain][loc]['main_min'],
                        dict_bd[chain][loc]['main_mean'],
                        dict_bd[chain][loc]['main_max'],

                        dict_bd[chain][loc]['side_min_surface'],
                        dict_bd[chain][loc]['side_mean_surface'],
                        dict_bd[chain][loc]['side_max_surface'],
                        dict_bd[chain][loc]['main_min_surface'],
                        dict_bd[chain][loc]['main_mean_surface'],
                        dict_bd[chain][loc]['main_max_surface'],
                        ])) + ','
                else:
                    str_output = str_output + ',,,,,,,,,,,,'

                #Burial Distance
                if self.dict_workflow['interface_distance']:
                    str_output = str_output + ','.join(map(str, [
                        dict_int_dist[chain][loc]['side_min'],                                   
                        dict_int_dist[chain][loc]['side_mean'],
                        dict_int_dist[chain][loc]['side_max'],
                        dict_int_dist[chain][loc]['main_min'],
                        dict_int_dist[chain][loc]['main_mean'],
                        dict_int_dist[chain][loc]['main_max'],
                        ])) + ','
                else:
                    str_output = str_output + ',,,,,,'

                #Newline
                str_output = str_output + '\n'            

            #Output our file
            with open(self.directory + self.output_prefix + '_structure.csv', 'w') as file_output:
                file_output.write(str_output)

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
