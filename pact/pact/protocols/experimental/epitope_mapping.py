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

"""epitope_mapping - method to indentify a epitope from a deep sequencing dataset"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Protocols:Epitope Mapping Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from configparser import ConfigParser, NoSectionError, NoOptionError
from pact.pact_common import file_checker
from sys import platform
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class epitope_mapping:
    """Calculate a putative epitope and compare epitope labellings"""

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
            print("[Protocols:Epitope Mapping Error] The config file is incorrect.")
            print("[Protocols:Epitope Mapping Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:Epitope Mapping Error] The config file is incorrect.")
            print("[Protocols:Epitope Mapping Error] There is something wrong with the [workflow] section (spelling???).")
            quit()

        #Get the PACT preamble
        self.pact_preamble = preamble

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Add trailing slash on directory
        if platform.startswith('win'):            
            if self.directory[-1:] != "\\":
                self.directory = self.directory + "\\"
        else:
            if self.directory[-1:] != "/":
                self.directory = self.directory + "/"
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
            Pact Combine Section (Required)
            *****************************************
            """
            #Check to see if the section is there
            if not self.obj_cfgparser.has_section('combinepact'):           
                print("[Protocols:Epitope Mapping Error] The combinepact config file is incorrect.")
                print("[Protocols:Epitope Mapping Error] There is something wrong with the [combinepact] section.")
                quit()

            #Import our combinepact class
            try:
                from pact.analysis.combine_pact import combine_pact
            except ImportError:
                print("[Protocols:Epitope Mapping Error] combine_pact was not found.")

            #Create the object then call the merger
            obj_combine = combine_pact(self.obj_cfgparser, self.dict_programs, {})

            #The dict will be like {'dataset name': {data...
            dict_merged_datasets = obj_combine.combine_pact()

            """
            *****************************************
            PDB Import Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['pdb_import']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('pdb_import'):           
                    print("[Protocols:Epitope Mapping Error] The pdb_import config file is incorrect.")
                    print("[Protocols:Epitope Mapping Error] There is something wrong with the [pdb_import] section.")
                    quit()

                #Import our combinepact class
                try:
                    from pact.analysis.pdb_import import pdb_import
                except ImportError:
                    print("[Protocols:Epitope Mapping Error] pdb_import was not found.")

                #Create the object then call the merger
                obj_pdb = pdb_import(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #The dict will be like {'pdb name': {data...
                dict_pdb = obj_pdb.pdb_import()

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:Epitope Mapping Error] This protocol needs to be ran within the context of PACT.")
