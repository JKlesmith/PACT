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

"""PACT vs Analysis - Compare two PACT datasets against an analysis"""

from sys import version_info

#Setup our protocol name and version
str_protocol_name = "PACT vs Analysis"

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

class pact_vs_analysis:
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
            T-Test of Two Groups
            *****************************************
            """
            if self.dict_workflow['aa_compare_ttest']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('aa_compare_ttest'):           
                    print("[Protocols:Analysis Error] The aa_compare_ttest config file is incorrect.")
                    print("[Protocols:Analysis Error] There is something wrong with the [aa_compare_ttest] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.sequence.aa_fitmet_compare import aa_fitmet_compare
                except ImportError:
                    print("[Protocols:Analysis Error] aa_fitmet_compare was not found.")
        
                #Create the object then call the merger
                obj_aac = aa_fitmet_compare(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #Run the main routine
                print("[Protocols:Analysis] T-Test of amino acid groups")
                file_output.write(obj_aac.aa_fitmet_compare(dict_merged_datasets))

            """
            *****************************************
            Count our mutations
            *****************************************
            """
            if self.dict_workflow['threshold_count']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('threshold_count'):           
                    print("[Protocols:Analysis Error] The threshold_count config file is incorrect.")
                    print("[Protocols:Analysis Error] There is something wrong with the [threshold_count] section.")
                    quit()

                #Create our object
                try:
                    from pact.analysis.sequence.threshold_count import threshold_count
                except ImportError:
                    print("[Protocols:Analysis Error] threshold_count was not found.")

                #Create the object then call the analysis
                obj_tc = threshold_count(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #Count them
                print("[Protocols:Analysis] Count of mutations above and below a cutoff")
                file_output.write(obj_tc.threshold_count(dict_merged_datasets))

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
