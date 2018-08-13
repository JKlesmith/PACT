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

"""PACT vs PACT - Compare two PACT datasets against each other"""

from sys import version_info

#Setup our protocol name and version
str_protocol_name = "PACT vs PACT"

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Protocols:" + str_protocol_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from configparser import ConfigParser, NoSectionError, NoOptionError
from os import access, W_OK, R_OK
from time import strftime

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class pact_vs_pact:
    """Provide an entrypoint for multiple analyses"""

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
            print("[Protocols:Analysis Error] The config file is incorrect.")
            print("[Protocols:Analysis Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:Analysis Error] The config file is incorrect.")
            print("[Protocols:Analysis Error] There is something wrong with the [workflow] section (spelling???).")
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
        """Provide a protocol that does general analyses that don't need a full protocol"""

        #Create a output log file that we can append to
        with open(self.directory + self.output_prefix + "_" + strftime("%m_%d_%Y") + "-" +
                 strftime("%H_%M_%S") + '_output.txt', 'w') as file_output:
            file_output.write(self.pact_preamble + "\n")

            """
            *****************************************
            Pact Combine (Required)
            *****************************************
            """
            #Check to see if the section is there
            if not self.obj_cfgparser.has_section('combinepact'):           
                print("[Protocols:Analysis Error] The combinepact config file is incorrect.")
                print("[Protocols:Analysis Error] There is something wrong with the [combinepact] section.")
                quit()

            #Import our combinepact class
            try:
                from pact.analysis.combine_pact import combine_pact
            except ImportError:
                print("[Protocols:Analysis Error] combine_pact was not found.")

            #Create the object then call the merger
            obj_combine = combine_pact(self.obj_cfgparser, self.dict_programs, {})

            #Print Section Progress
            print("[Protocols:Analysis] Combine PACT")

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
                    print("[Protocols:Analysis Error] The pdb_import config file is incorrect.")
                    print("[Protocols:Analysis Error] There is something wrong with the [pdb_import] section.")
                    quit()

                #Import our combinepact class
                try:
                    from pact.analysis.pdb_import import pdb_import
                except ImportError:
                    print("[Protocols:Analysis Error] pdb_import was not found.")

                #Create the object then call the merger
                obj_pdb = pdb_import(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #Print Section Progress
                print("[Protocols:Analysis] PDB Import")

                #The dict will be like {'pdb name': {data...
                dict_pdb = obj_pdb.pdb_import()

            """
            *****************************************
            Assign colors to classifiers
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['classifier_color']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('classifier_color'):           
                    print("[Protocols:Analysis Error] The classifier_color config file is incorrect.")
                    print("[Protocols:Analysis Error] There is something wrong with the [classifier_color] section.")
                    quit()

                #Import our combinepact class
                try:
                    from pact.analysis.classifier_color import classifier_color
                except ImportError:
                    print("[Protocols:Analysis Error] classifer_color was not found.")

                #Create the object then call the merger
                obj_classcolor = classifier_color(self.obj_cfgparser, self.dict_programs, {})

                #Print Section Progress
                print("[Protocols:Analysis] Classifier Color")

                #This section returns a dict of [loc][mut] = "color"
                if self.obj_cfgparser.get("classifier_color", "classifier").split(',')[0] == "pdb":
                    dict_custom_color = obj_classcolor.classifier_color(dict_merged_datasets, dict_pdb, "pdb")

            """
            *****************************************
            Set vs Set Section
            *****************************************
            """
            #Check to see if the section is there
            if not self.obj_cfgparser.has_section('setvsset'):           
                print("[Protocols:Analysis Error] The setvsset config file is incorrect.")
                print("[Protocols:Analysis Error] There is something wrong with the [setvsset] section.")
                quit()

            #Import our setvsset class
            try:
                from pact.analysis.set_vs_set import set_vs_set
            except ImportError:
                print("[Protocols:Analysis Error] set_vs_set was not found.")

            #Create the object then call the merger
            obj_svs = set_vs_set(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

            #Do we have structural data?
            if self.dict_workflow['classifier_color']:
                print("[Protocols:Analysis] Dataset vs Dataset")
                file_output.write(obj_svs.set_vs_set(dict_merged_datasets, dict_custom_color))
            else:
                print("[Protocols:Analysis] Dataset vs Dataset")
                file_output.write(obj_svs.set_vs_set(dict_merged_datasets))

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
