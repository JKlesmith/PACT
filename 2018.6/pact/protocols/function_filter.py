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

"""enzyme solubility filter - protocol to filter mutations that pass our enzyme solubility criteria"""

from sys import version_info

#Setup our protocol name
str_protocol_name = "Enzyme Filter"

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
from pact.pact_common import file_checker, pretty_counter_dicts, open_pact_file
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class enzyme_filter:
    """Calculate mutations on a enzyme solubility filter"""

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

    def classifier_count(self, dict_classifiers):
        """Count our classifiers"""

        """
        *****************************************
        Import Classifiers
        *****************************************
        """
        #Only import and run if selected
        if self.dict_workflow['bayes_count']:

            #Check to see if the section is there
            if not self.obj_cfgparser.has_section('bayes_count'):           
                print("[Protocols:Enzyme Solubility Filter Error] The bayes_count config file is incorrect.")
                print("[Protocols:Enzyme Solubility Filter Error] There is something wrong with the [bayes_count] section.")
                quit()

            #Import our combinepact class
            try:
                from pact.classification.bayes_count import bayes_count
            except ImportError:
                print("[Protocols:Enzyme Solubility Filter Error] bayes_count was not found.")

            #Create the object then call the merger
            obj_bayes = bayes_count(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

            #Print Section Progress
            print("[Protocols:Enzyme Solubility Filter] Bayes Count")

            print("Basal Rate")
            list_ben, list_neu, list_del = obj_bayes.bayes_count(dict_classifiers)

        """
        *****************************************
        Count Section
        *****************************************
        """


        """
        print("")
        #PSSM
        print("PSSM >=3")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_pssm', 3, ">=", 20, "<")

        print("PSSM < 3 & >= 0")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_pssm', 0, ">=", 2, "<")

        print("PSSM < 0")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_pssm', -20, ">=", 0, "<")

        print("")
            
        #Fraction Burial
        print("frac_burial 0 < 30")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'frac_burial', 0, ">=", 0.3, "<")

        print("frac_burial 30 < 60")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'frac_burial', 0.3, ">=", 0.6, "<")

        print("frac_burial 60 < 90")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'frac_burial', 0.6, ">=", 0.9, "<")

        print("frac_burial 90 to 100")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'frac_burial', 0.9, ">=", 1, "<=")

        print("")
            
        #Contact Number
        print("contact_number 0 < 10")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'contact_number', 0, ">=", 10, "<")

        print("contact_number 10 < 20")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'contact_number', 10, ">=", 20, "<")

        print("contact_number 20+")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'contact_number', 20, ">=", 200, "<")

        print("")
            
        #Dist to Active
        print("dist_to_active 0 < 5")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 0, ">=", 5, "<")

        print("dist_to_active 5 < 10")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 5, ">=", 10, "<")

        print("dist_to_active 10 < 15")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 10, ">=", 15, "<")

        print("dist_to_active 15 < 20")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 15, ">=", 20, "<")

        print("dist_to_active 20 < 25")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 20, ">=", 25, "<")

        print("dist_to_active 25 < 30")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 25, ">=", 30, "<")

        print("dist_to_active 30+")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'dist_to_active', 30, ">=", 200, "<")

        print("")
            
        #Wild-Type Percentage at Residue
        print("wt_percent 0 to 15")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_percent', 0, ">=", 16, "<")

        print("wt_percent 16 to 30")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_percent', 16, ">=", 31, "<")

        print("wt_percent 31 to 45")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_percent', 31, ">=", 45, "<")

        print("wt_percent 46 to 60")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_percent', 45, ">=", 60, "<")

        print("wt_percent 60 to 100")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_percent', 60, ">=", 101, "<")

        print("")
        #Max Site Percentage
        print("max_percent 0 to 15")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'max_percent', 0, ">=", 16, "<")

        print("max_percent 16 to 30")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'max_percent', 16, ">=", 31, "<")

        print("max_percent 31 to 45")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'max_percent', 31, ">=", 45, "<")

        print("max_percent 46 to 60")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'max_percent', 45, ">=", 60, "<")

        print("max_percent 60 to 100")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'max_percent', 60, ">=", 101, "<")

        print("")
        #Mutation Percentage
        print("mut_percent 0")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_percent', 0, "=", 0, "=")

        print("mut_percent 1 to 15")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_percent', 1, ">=", 16, "<")

        print("mut_percent 16 to 30")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_percent', 16, ">=", 31, "<")

        print("mut_percent 31 to 45")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_percent', 31, ">=", 45, "<")

        print("mut_percent 46 to 60")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_percent', 45, ">=", 60, "<")

        print("mut_percent 60 to 100")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'mut_percent', 60, ">=", 101, "<")

        print("")
        #WT Consensus?
        print("wt_max_percent 0")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_max_percent', 0, "=", 0, "=")

        print("wt_max_percent 1")
        obj_bayes.classifier_count_range(list_ben, list_neu, list_del, 'wt_max_percent', 1, "=", 1, "=")

        print("")
        #Mutation Consensus?
        print("mut cons? 0")
        obj_bayes.classifier_count_keycompare(list_ben, list_neu, list_del, 'mut_percent', 'max_percent', "!")

        print("mut cons? 1")
        obj_bayes.classifier_count_keycompare(list_ben, list_neu, list_del, 'mut_percent', 'max_percent', "=")
            

        print("")
        #Size Change
        print("Size: SB, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.85, '>=', 1, '<=', ['SB'])

        print("Size: BS, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.85, '>=', 1, '<=', ['BS'])

        print("Size: SS/BB, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.85, '>=', 1, '<=', ['BB', 'SS'])

        print("Size: PX/XP, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.85, '>=', 1, '<=', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: SB, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.5, '>=', 0.85, '<', ['SB'])

        print("Size: BS, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.5, '>=', 0.85, '<', ['BS'])

        print("Size: SS/BB, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.5, '>=', 0.85, '<', ['BB', 'SS'])

        print("Size: PX/XP, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0.5, '>=', 0.85, '<', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: SB, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0, '>=', 0.5, '<', ['SB'])

        print("Size: BS, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0, '>=', 0.5, '<', ['BS'])

        print("Size: SS/BB, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0, '>=', 0.5, '<', ['BB', 'SS'])

        print("Size: PX/XP, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0, '>=', 0.5, '<', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: ALL")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'size', 
                                                0, '>=', 1, '<=', ['SP', 'BP', 'PS', 'PB'])


        print("Size: Contact Number")
        #Size Change
        print("Size: SB, Contact Number >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                20, '>=', 100, '<=', ['SB'])

        print("Size: BS, Contact Number >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                20, '>=', 100, '<=', ['BS'])

        print("Size: SS/BB, Contact Number >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                20, '>=', 100, '<=', ['BB', 'SS'])

        print("Size: PX/XP, Contact Number >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                20, '>=', 100, '<=', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: SB, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                11, '>=', 20, '<', ['SB'])

        print("Size: BS, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                11, '>=', 20, '<', ['BS'])

        print("Size: SS/BB, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                11, '>=', 20, '<', ['BB', 'SS'])

        print("Size: PX/XP, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                11, '>=', 20, '<', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: SB, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                0, '>=', 11, '<', ['SB'])

        print("Size: BS, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                0, '>=', 11, '<', ['BS'])

        print("Size: SS/BB, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                0, '>=', 11, '<', ['BB', 'SS'])

        print("Size: PX/XP, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                0, '>=', 11, '<', ['SP', 'BP', 'PS', 'PB'])

            
        print("Size: Distance to Active")
        #Size Change
        print("Size: SB, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                0, '>=', 10, '<', ['SB'])

        print("Size: BS, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                0, '>=', 10, '<', ['BS'])

        print("Size: SS/BB, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                0, '>=', 10, '<', ['BB', 'SS'])

        print("Size: PX/XP, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                0, '>=', 10, '<', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: SB, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                10, '>=', 20, '<', ['SB'])

        print("Size: BS, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                10, '>=', 20, '<', ['BS'])

        print("Size: SS/BB, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                10, '>=', 20, '<', ['BB', 'SS'])

        print("Size: PX/XP, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                10, '>=', 20, '<', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Size: SB, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                20, '>=', 100, '<', ['SB'])

        print("Size: BS, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                20, '>=', 100, '<', ['BS'])

        print("Size: SS/BB, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                20, '>=', 100, '<', ['BB', 'SS'])

        print("Size: PX/XP, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'size', 
                                                20, '>=', 100, '<', ['SP', 'BP', 'PS', 'PB'])

            

        print("")
        print("HP/HL: Frac Burial")
        #Size Change
        print("HP/HL: No Change, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0.85, '>=', 1, '<=', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0.85, '>=', 1, '<=', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Frac Burial >= 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0.85, '>=', 1, '<=', ['LB', 'LC'])

        print("")
        print("HP/HL: No Change, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0.5, '>=', 0.85, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0.5, '>=', 0.85, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Frac Burial >=0.5 < 0.85")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0.5, '>=', 0.85, '<', ['LB', 'LC'])

        print("")
        print("HP/HL: No Change, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0, '>=', 0.5, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0, '>=', 0.5, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Frac Burial < 0.5")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'frac_burial', 'philic_phobic', 
                                                0, '>=', 0.5, '<', ['LB', 'LC'])


        print("")
        print("HP/HL: Contact Number")
        #Size Change
        print("HP/HL: No Change, Frac Burial >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                20, '>=', 100, '<=', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Frac Burial >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                20, '>=', 100, '<=', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Frac Burial >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                20, '>=', 100, '<=', ['LB', 'LC'])

        print("")
        print("HP/HL: No Change, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                11, '>=', 20, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                11, '>=', 20, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Frac Burial >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                11, '>=', 20, '<', ['LB', 'LC'])

        print("")
        print("HP/HL: No Change, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                0, '>=', 11, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                0, '>=', 11, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Frac Burial <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'philic_phobic', 
                                                0, '>=', 11, '<', ['LB', 'LC'])


        print("")
        print("HP/HL: Distance to Active")
        #Size Change
        print("HP/HL: No Change, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                0, '>=', 10, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                0, '>=', 10, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Distance to Active <10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                0, '>=', 10, '<', ['LB', 'LC'])


        print("")
        print("HP/HL: No Change, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                10, '>=', 20, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                10, '>=', 20, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Distance to Active 10 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                10, '>=', 20, '<', ['LB', 'LC'])

        print("")
        print("HP/HL: No Change, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                20, '>=', 100, '<', ['BB', 'LL', 'CB', 'BC', 'CC'])

        print("HP/HL: Hydrophobic to Hydrophilic, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                20, '>=', 100, '<', ['BL', 'CL'])

        print("HP/HL: Hydrophilic to Hydrophobic, Distance to Active >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'dist_to_active', 'philic_phobic', 
                                                20, '>=', 100, '<', ['LB', 'LC'])
        """






        print("Contact Number: Non-proline vs proline")
        #Size Change
        print("Non-Pro Contact Number >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                20, '>=', 100, '<=', ['SB', 'BS', 'BB', 'SS'])

        print("PRO Contact Number >= 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                20, '>=', 100, '<=', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Non-Pro Contact Number >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                11, '>=', 20, '<', ['SB', 'BS', 'BB', 'SS'])

        print("PRO Contact Number >=11 < 20")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                11, '>=', 20, '<', ['SP', 'BP', 'PS', 'PB'])

        print("")
        print("Non-Pro Contact Number <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                0, '>=', 11, '<', ['SB', 'BS', 'BB', 'SS'])

        print("PRO Contact Number <= 10")
        obj_bayes.classifier_count_rangelist(list_ben, list_neu, list_del, 'contact_number', 'size', 
                                                0, '>=', 11, '<', ['SP', 'BP', 'PS', 'PB'])


        return

    def bayes_model_score(self, dict_classifiers):
        """Score our model with a dataset"""

        #Setup our prob(BEN, NEU, DEL) with our priors
        p_ben = 0.02280
        p_neu = 0.18678
        p_del = 0.79042

        #Setup our evidence with 1 (since we're multiplying)
        p_evidence = 1

        if self.obj_cfgparser.get('bayes_model_score', 'pssm_variant').lower() == "true":
            #Setup the PSSM
            if int(dict_classifiers['mut_pssm']) >= 3:
                p_ben *= 0.24306
                p_neu *= 0.10987
                p_del *= 0.02953
                p_evidence *= 0.04615
            elif int(dict_classifiers['mut_pssm']) >= 0 and int(dict_classifiers['mut_pssm']) < 3:
                p_ben *= 0.36806
                p_neu *= 0.34395
                p_del *= 0.12511
                p_evidence *= 0.16166
            else:
                p_ben *= 0.38889
                p_neu *= 0.54618
                p_del *= 0.84536
                p_evidence *= 0.74618

        if self.obj_cfgparser.get('bayes_model_score', 'frac_burial').lower() == "true":
            #Setup the fraction burial
            if float(dict_classifiers['frac_burial']) < 0.3:
                p_ben *= 0.05988
                p_neu *= 0.15716
                p_del *= 0.02487
                p_evidence *= 0.05038
            elif float(dict_classifiers['frac_burial']) >= 0.3 and float(dict_classifiers['frac_burial']) < 0.6:
                p_ben *= 0.28743
                p_neu *= 0.28289
                p_del *= 0.07186
                p_evidence *= 0.11619
            elif float(dict_classifiers['frac_burial']) >= 0.6 and float(dict_classifiers['frac_burial']) < 0.9:
                p_ben *= 0.35928
                p_neu *= 0.31579
                p_del *= 0.20176
                p_evidence *= 0.22665
            else:
                p_ben *= 0.29341
                p_neu *= 0.24415
                p_del *= 0.70150
                p_evidence *= 0.60677
       
        if self.obj_cfgparser.get('bayes_model_score', 'contact_number').lower() == "true":
            #Setup the contact number
            if int(dict_classifiers['contact_number']) < 10:
                p_ben *= 0.13174
                p_neu *= 0.13670
                p_del *= 0.02574
                p_evidence *= 0.04888
            elif int(dict_classifiers['contact_number']) >= 11 and int(dict_classifiers['contact_number']) < 20:
                p_ben *= 0.61078
                p_neu *= 0.68202
                p_del *= 0.39074
                p_evidence *= 0.45016
            else:
                p_ben *= 0.25749
                p_neu *= 0.18129
                p_del *= 0.58352
                p_evidence *= 0.50096

        if self.obj_cfgparser.get('bayes_model_score', 'wt_cons').lower() == "true":
            #Setup wt consen?
            if int(dict_classifiers['wt_max_percent']) == 1:
                p_ben *= 0.25749
                p_neu *= 0.51462
                p_del *= 0.70599
                p_evidence *= 0.66002
            else:
                p_ben *= 0.74251
                p_neu *= 0.48538
                p_del *= 0.29401
                p_evidence *= 0.33998

        if self.obj_cfgparser.get('bayes_model_score', 'variant_cons').lower() == "true":
            #Setup variant consen?
            if int(dict_classifiers['mut_percent']) == int(dict_classifiers['max_percent']):
                p_ben *= 0.14970
                p_neu *= 0.04459
                p_del *= 0.00967
                p_evidence *= 0.01939
            else:
                p_ben *= 0.85030
                p_neu *= 0.95541
                p_del *= 0.99033
                p_evidence *= 0.98061

        if self.obj_cfgparser.get('bayes_model_score', 'd_to_a').lower() == "true":
            #Setup the distance to active site
            if float(dict_classifiers['dist_to_active']) < 5:
                p_ben *= 0.01796
                p_neu *= 0.02705
                p_del *= 0.04713
                p_evidence *= 0.04014
            elif float(dict_classifiers['dist_to_active']) >= 5 and float(dict_classifiers['dist_to_active']) < 10:
                p_ben *= 0.10778
                p_neu *= 0.06287
                p_del *= 0.17796
                p_evidence *= 0.14514
            elif float(dict_classifiers['dist_to_active']) >= 10 and float(dict_classifiers['dist_to_active']) < 15:
                p_ben *= 0.11976
                p_neu *= 0.14474
                p_del *= 0.21470
                p_evidence *= 0.18774
            elif float(dict_classifiers['dist_to_active']) >= 15 and float(dict_classifiers['dist_to_active']) < 20:
                p_ben *= 0.24551
                p_neu *= 0.23684
                p_del *= 0.26851
                p_evidence *= 0.24741
            elif float(dict_classifiers['dist_to_active']) >= 20 and float(dict_classifiers['dist_to_active']) < 25:
                p_ben *= 0.20359
                p_neu *= 0.27924
                p_del *= 0.19503
                p_evidence *= 0.20030
            elif float(dict_classifiers['dist_to_active']) >= 25 and float(dict_classifiers['dist_to_active']) < 30:
                p_ben *= 0.24551
                p_neu *= 0.21491
                p_del *= 0.08035
                p_evidence *= 0.10486
            else:
                p_ben *= 0.05988
                p_neu *= 0.03436
                p_del *= 0.01633
                p_evidence *= 0.01980

        if self.obj_cfgparser.get('bayes_model_score', 'mut_percent').lower() == "true":
            #Setup mutation percentage
            if int(dict_classifiers['mut_percent']) == 0:
                p_ben *= 0.19760
                p_neu *= 0.34649
                p_del *= 0.64605
                p_evidence *= 0.57987
            elif int(dict_classifiers['mut_percent']) >= 1 and int(dict_classifiers['mut_percent']) < 15:
                p_ben *= 0.54491
                p_neu *= 0.54605
                p_del *= 0.32614
                p_evidence *= 0.37220
            elif int(dict_classifiers['mut_percent']) >= 15 and int(dict_classifiers['mut_percent']) < 30:
                p_ben *= 0.14371
                p_neu *= 0.08114
                p_del *= 0.02038
                p_evidence *= 0.03454
            elif int(dict_classifiers['mut_percent']) >= 30 and int(dict_classifiers['mut_percent']) < 45:
                p_ben *= 0.05988
                p_neu *= 0.01681
                p_del *= 0.00415
                p_evidence *= 0.00778
            elif int(dict_classifiers['mut_percent']) >= 45 and int(dict_classifiers['mut_percent']) < 60:
                p_ben *= 0.03593
                p_neu *= 0.00585
                p_del *= 0.00173
                p_evidence *= 0.00328
            else:
                p_ben *= 0.01796
                p_neu *= 0.00365
                p_del *= 0.00155
                p_evidence *= 0.00232      

        if self.obj_cfgparser.get('bayes_model_score', 'max_percent').lower() == "true":
            #Setup wt %
            if int(dict_classifiers['max_percent']) < 15:
                p_ben *= 0.01796
                p_neu *= 0.01023
                p_del *= 0.00535
                p_evidence *= 0.00655
            elif int(dict_classifiers['max_percent']) >= 15 and int(dict_classifiers['max_percent']) < 30:
                p_ben *= 0.49701
                p_neu *= 0.47076
                p_del *= 0.16721
                p_evidence *= 0.23143
            elif int(dict_classifiers['max_percent']) >= 30 and int(dict_classifiers['max_percent']) < 45:
                p_ben *= 0.20359
                p_neu *= 0.20833
                p_del *= 0.22456
                p_evidence *= 0.22105
            elif int(dict_classifiers['max_percent']) >= 45 and int(dict_classifiers['max_percent']) < 60:
                p_ben *= 0.18563
                p_neu *= 0.12135
                p_del *= 0.17775
                p_evidence *= 0.16739
            else:
                p_ben *= 0.09581
                p_neu *= 0.18933
                p_del *= 0.42512
                p_evidence *= 0.37357
            
        if self.obj_cfgparser.get('bayes_model_score', 'wt_percent').lower() == "true":
            #Setup wt %
            if int(dict_classifiers['wt_percent']) < 15:
                p_ben *= 0.59880
                p_neu *= 0.30482
                p_del *= 0.13750
                p_evidence *= 0.17927
            elif int(dict_classifiers['wt_percent']) >= 15 and int(dict_classifiers['wt_percent']) < 30:
                p_ben *= 0.25150
                p_neu *= 0.36623
                p_del *= 0.20332
                p_evidence *= 0.23484
            elif int(dict_classifiers['wt_percent']) >= 30 and int(dict_classifiers['wt_percent']) < 45:
                p_ben *= 0.03593
                p_neu *= 0.09576
                p_del *= 0.14320
                p_evidence *= 0.13190
            elif int(dict_classifiers['wt_percent']) >= 45 and int(dict_classifiers['wt_percent']) < 60:
                p_ben *= 0.04192
                p_neu *= 0.07310
                p_del *= 0.13336
                p_evidence *= 0.12002
            else:
                p_ben *= 0.07186
                p_neu *= 0.16009
                p_del *= 0.38262
                p_evidence *= 0.33397

        if self.obj_cfgparser.get('bayes_model_score', 'pro_v_contactnum').lower() == "true":
            #PRO/NON-PRO vs Contact Number
            #Non-Pro
            if int(dict_classifiers['contact_number']) >= 20 and dict_classifiers['size'] in ['SB', 'BS', 'BB', 'SS']:
                p_ben *= 0.9940
                p_neu *= 0.9435
                p_del *= 0.9156
                p_evidence *= 0.4771
            #Pro
            elif int(dict_classifiers['contact_number']) >= 20 and dict_classifiers['size'] in ['SP', 'BP', 'PS', 'PB']:
                p_ben *= 0.0060
                p_neu *= 0.0565
                p_del *= 0.0844
                p_evidence *= 0.0410
            #Non-Pro
            elif (int(dict_classifiers['contact_number']) >= 11 and 
                int(dict_classifiers['contact_number']) < 20 and 
                dict_classifiers['size'] in ['SB', 'BS', 'BB', 'SS']):
                p_ben *= 0.9485
                p_neu *= 0.9708
                p_del *= 0.8719
                p_evidence *= 0.3835
            #Pro
            elif (int(dict_classifiers['contact_number']) >= 11 and
                 int(dict_classifiers['contact_number']) < 20 and 
                 dict_classifiers['size'] in ['SP', 'BP', 'PS', 'PB']):
                p_ben *= 0.0515
                p_neu *= 0.0292
                p_del *= 0.1281
                p_evidence *= 0.0425
            #Non-Pro
            elif (int(dict_classifiers['contact_number']) < 11 and 
                dict_classifiers['size'] in ['SB', 'BS', 'BB', 'SS']):
                p_ben *= 0.7778
                p_neu *= 0.9231
                p_del *= 0.8373
                p_evidence *= 0.0644
            #Pro
            elif (int(dict_classifiers['contact_number']) < 11 and 
                 dict_classifiers['size'] in ['SP', 'BP', 'PS', 'PB']):
                p_ben *= 0.2222
                p_neu *= 0.0769
                p_del *= 0.1627
                p_evidence *= 0.0086
        
        #Calculate our probs
        prob_ben = p_ben / p_evidence
        prob_neu = p_neu / p_evidence
        prob_del = p_del / p_evidence
        
        return prob_ben, prob_neu, prob_del



    def strict_filter_old(self, dict_classifiers):

        if (int(dict_classifiers['mut_pssm']) >= 0 and 
            float(dict_classifiers['dist_to_active']) >= 15 and
            int(dict_classifiers['contact_number']) <= 16 and
            dict_classifiers['size'] not in ['SP', 'BP', 'PS', 'PB']):
            mut_pass_fail = "BEN"
        else:
            mut_pass_fail = "DEL"

        return mut_pass_fail

    def strict_filter_new(self, dict_classifiers):


        #int(dict_classifiers['wt_max_percent']) != 1 and
        #int(dict_classifiers['mut_pssm']) >= 0 and
        #int(dict_classifiers['wt_max_percent']) < int(dict_classifiers['mut_percent']) and

        if (int(dict_classifiers['mut_pssm']) >= 0 and
            float(dict_classifiers['dist_to_active']) >= 15 and
            int(dict_classifiers['contact_number']) <= 16 and
            float(dict_classifiers['frac_burial']) < 0.4 and
            int(dict_classifiers['wt_max_percent']) != 1 and
            dict_classifiers['size'] not in ['SP', 'BP', 'PS', 'PB']):
            mut_pass_fail = "BEN"
        else:
            mut_pass_fail = "DEL"

        return mut_pass_fail



    def bayes_model_score_combo(self, dict_classifiers, dict_options):
        """Score our model with a dataset"""

        #Setup our prob(BEN, NEU, DEL) with our priors
        p_ben = 0.02280
        p_neu = 0.18678
        p_del = 0.79042

        #Setup our evidence with 1 (since we're multiplying)
        p_evidence = 1

        if dict_options[0]:
            #Setup the PSSM
            if int(dict_classifiers['mut_pssm']) >= 3:
                p_ben *= 0.24306
                p_neu *= 0.10987
                p_del *= 0.02953
                p_evidence *= 0.04615
            elif int(dict_classifiers['mut_pssm']) >= 0 and int(dict_classifiers['mut_pssm']) < 3:
                p_ben *= 0.36806
                p_neu *= 0.34395
                p_del *= 0.12511
                p_evidence *= 0.16166
            else:
                p_ben *= 0.38889
                p_neu *= 0.54618
                p_del *= 0.84536
                p_evidence *= 0.74618

        if dict_options[1]:
            #Setup the fraction burial
            if float(dict_classifiers['frac_burial']) < 0.3:
                p_ben *= 0.05988
                p_neu *= 0.15716
                p_del *= 0.02487
                p_evidence *= 0.05038
            elif float(dict_classifiers['frac_burial']) >= 0.3 and float(dict_classifiers['frac_burial']) < 0.6:
                p_ben *= 0.28743
                p_neu *= 0.28289
                p_del *= 0.07186
                p_evidence *= 0.11619
            elif float(dict_classifiers['frac_burial']) >= 0.6 and float(dict_classifiers['frac_burial']) < 0.9:
                p_ben *= 0.35928
                p_neu *= 0.31579
                p_del *= 0.20176
                p_evidence *= 0.22665
            else:
                p_ben *= 0.29341
                p_neu *= 0.24415
                p_del *= 0.70150
                p_evidence *= 0.60677
       
        if dict_options[2]:
            #Setup the contact number
            if int(dict_classifiers['contact_number']) < 10:
                p_ben *= 0.13174
                p_neu *= 0.13670
                p_del *= 0.02574
                p_evidence *= 0.04888
            elif int(dict_classifiers['contact_number']) >= 11 and int(dict_classifiers['contact_number']) < 20:
                p_ben *= 0.61078
                p_neu *= 0.68202
                p_del *= 0.39074
                p_evidence *= 0.45016
            else:
                p_ben *= 0.25749
                p_neu *= 0.18129
                p_del *= 0.58352
                p_evidence *= 0.50096

        if dict_options[3]:
            #Setup wt consen?
            if int(dict_classifiers['wt_max_percent']) == 1:
                p_ben *= 0.25749
                p_neu *= 0.51462
                p_del *= 0.70599
                p_evidence *= 0.66002
            else:
                p_ben *= 0.74251
                p_neu *= 0.48538
                p_del *= 0.29401
                p_evidence *= 0.33998

        if dict_options[4]:
            #Setup variant consen?
            if int(dict_classifiers['mut_percent']) == int(dict_classifiers['max_percent']):
                p_ben *= 0.14970
                p_neu *= 0.04459
                p_del *= 0.00967
                p_evidence *= 0.01939
            else:
                p_ben *= 0.85030
                p_neu *= 0.95541
                p_del *= 0.99033
                p_evidence *= 0.98061

        if dict_options[5]:
            #Setup the distance to active site
            if float(dict_classifiers['dist_to_active']) < 5:
                p_ben *= 0.01796
                p_neu *= 0.02705
                p_del *= 0.04713
                p_evidence *= 0.04014
            elif float(dict_classifiers['dist_to_active']) >= 5 and float(dict_classifiers['dist_to_active']) < 10:
                p_ben *= 0.10778
                p_neu *= 0.06287
                p_del *= 0.17796
                p_evidence *= 0.14514
            elif float(dict_classifiers['dist_to_active']) >= 10 and float(dict_classifiers['dist_to_active']) < 15:
                p_ben *= 0.11976
                p_neu *= 0.14474
                p_del *= 0.21470
                p_evidence *= 0.18774
            elif float(dict_classifiers['dist_to_active']) >= 15 and float(dict_classifiers['dist_to_active']) < 20:
                p_ben *= 0.24551
                p_neu *= 0.23684
                p_del *= 0.26851
                p_evidence *= 0.24741
            elif float(dict_classifiers['dist_to_active']) >= 20 and float(dict_classifiers['dist_to_active']) < 25:
                p_ben *= 0.20359
                p_neu *= 0.27924
                p_del *= 0.19503
                p_evidence *= 0.20030
            elif float(dict_classifiers['dist_to_active']) >= 25 and float(dict_classifiers['dist_to_active']) < 30:
                p_ben *= 0.24551
                p_neu *= 0.21491
                p_del *= 0.08035
                p_evidence *= 0.10486
            else:
                p_ben *= 0.05988
                p_neu *= 0.03436
                p_del *= 0.01633
                p_evidence *= 0.01980

        if dict_options[6]:
            #Setup mutation percentage
            if int(dict_classifiers['mut_percent']) == 0:
                p_ben *= 0.19760
                p_neu *= 0.34649
                p_del *= 0.64605
                p_evidence *= 0.57987
            elif int(dict_classifiers['mut_percent']) >= 1 and int(dict_classifiers['mut_percent']) < 15:
                p_ben *= 0.54491
                p_neu *= 0.54605
                p_del *= 0.32614
                p_evidence *= 0.37220
            elif int(dict_classifiers['mut_percent']) >= 15 and int(dict_classifiers['mut_percent']) < 30:
                p_ben *= 0.14371
                p_neu *= 0.08114
                p_del *= 0.02038
                p_evidence *= 0.03454
            elif int(dict_classifiers['mut_percent']) >= 30 and int(dict_classifiers['mut_percent']) < 45:
                p_ben *= 0.05988
                p_neu *= 0.01681
                p_del *= 0.00415
                p_evidence *= 0.00778
            elif int(dict_classifiers['mut_percent']) >= 45 and int(dict_classifiers['mut_percent']) < 60:
                p_ben *= 0.03593
                p_neu *= 0.00585
                p_del *= 0.00173
                p_evidence *= 0.00328
            else:
                p_ben *= 0.01796
                p_neu *= 0.00365
                p_del *= 0.00155
                p_evidence *= 0.00232      

        if dict_options[7]:
            #Setup max site percentage
            if int(dict_classifiers['max_percent']) < 15:
                p_ben *= 0.01796
                p_neu *= 0.01023
                p_del *= 0.00535
                p_evidence *= 0.00655
            elif int(dict_classifiers['max_percent']) >= 15 and int(dict_classifiers['max_percent']) < 30:
                p_ben *= 0.49701
                p_neu *= 0.47076
                p_del *= 0.16721
                p_evidence *= 0.23143
            elif int(dict_classifiers['max_percent']) >= 30 and int(dict_classifiers['max_percent']) < 45:
                p_ben *= 0.20359
                p_neu *= 0.20833
                p_del *= 0.22456
                p_evidence *= 0.22105
            elif int(dict_classifiers['max_percent']) >= 45 and int(dict_classifiers['max_percent']) < 60:
                p_ben *= 0.18563
                p_neu *= 0.12135
                p_del *= 0.17775
                p_evidence *= 0.16739
            else:
                p_ben *= 0.09581
                p_neu *= 0.18933
                p_del *= 0.42512
                p_evidence *= 0.37357
        
        if dict_options[8]:
            #Setup wt %
            if int(dict_classifiers['wt_percent']) < 15:
                p_ben *= 0.59880
                p_neu *= 0.30482
                p_del *= 0.13750
                p_evidence *= 0.17927
            elif int(dict_classifiers['wt_percent']) >= 15 and int(dict_classifiers['wt_percent']) < 30:
                p_ben *= 0.25150
                p_neu *= 0.36623
                p_del *= 0.20332
                p_evidence *= 0.23484
            elif int(dict_classifiers['wt_percent']) >= 30 and int(dict_classifiers['wt_percent']) < 45:
                p_ben *= 0.03593
                p_neu *= 0.09576
                p_del *= 0.14320
                p_evidence *= 0.13190
            elif int(dict_classifiers['wt_percent']) >= 45 and int(dict_classifiers['wt_percent']) < 60:
                p_ben *= 0.04192
                p_neu *= 0.07310
                p_del *= 0.13336
                p_evidence *= 0.12002
            else:
                p_ben *= 0.07186
                p_neu *= 0.16009
                p_del *= 0.38262
                p_evidence *= 0.33397
        
        if dict_options[9]:
            #PRO/NON-PRO vs Contact Number
            #Non-Pro
            if int(dict_classifiers['contact_number']) >= 20 and dict_classifiers['size'] in ['SB', 'BS', 'BB', 'SS']:
                p_ben *= 0.9940
                p_neu *= 0.9435
                p_del *= 0.9156
                p_evidence *= 0.4771
            #Pro
            elif int(dict_classifiers['contact_number']) >= 20 and dict_classifiers['size'] in ['SP', 'BP', 'PS', 'PB']:
                p_ben *= 0.0060
                p_neu *= 0.0565
                p_del *= 0.0844
                p_evidence *= 0.0410
            #Non-Pro
            elif (int(dict_classifiers['contact_number']) >= 11 and 
                int(dict_classifiers['contact_number']) < 20 and 
                dict_classifiers['size'] in ['SB', 'BS', 'BB', 'SS']):
                p_ben *= 0.9485
                p_neu *= 0.9708
                p_del *= 0.8719
                p_evidence *= 0.3835
            #Pro
            elif (int(dict_classifiers['contact_number']) >= 11 and
                 int(dict_classifiers['contact_number']) < 20 and 
                 dict_classifiers['size'] in ['SP', 'BP', 'PS', 'PB']):
                p_ben *= 0.0515
                p_neu *= 0.0292
                p_del *= 0.1281
                p_evidence *= 0.0425
            #Non-Pro
            elif (int(dict_classifiers['contact_number']) < 11 and 
                dict_classifiers['size'] in ['SB', 'BS', 'BB', 'SS']):
                p_ben *= 0.7778
                p_neu *= 0.9231
                p_del *= 0.8373
                p_evidence *= 0.0644
            #Pro
            elif (int(dict_classifiers['contact_number']) < 11 and 
                 dict_classifiers['size'] in ['SP', 'BP', 'PS', 'PB']):
                p_ben *= 0.2222
                p_neu *= 0.0769
                p_del *= 0.1627
                p_evidence *= 0.0086
        
        #Calculate our probs
        prob_ben = p_ben / p_evidence
        prob_neu = p_neu / p_evidence
        prob_del = p_del / p_evidence
        
        return prob_ben, prob_neu, prob_del

    def bayes_combo(self, dict_classifiers):

        from itertools import product
        bayes_combos = list(product([True,False], repeat=10))

        #Actual_Model
        output = ("BEN_BEN, BEN_NEU, BEN_DEL, NEU_BEN, NEU_NEU, NEU_DEL, DEL_BEN, DEL_NEU, DEL_DEL," +
                    "PSSM, Frac_Burial, CN, WT_Cons, Var_Cons, DtA, Mut_Perc, Max_Site_Perc, WT_Perc, PRO_CN\n")
        
        for combo in bayes_combos:
            #Setup our counters to measure correct, fp, fn
            list_count_ben = []
            list_count_neu = []
            list_count_del = []

            #Define our two datasets
            ds_activity = "amie_iso2_classified"
            #ds_activity = "LGK_Triple_classified"
            #ds_stability = "LGK_wt_classified"

            for loc in dict_classifiers:

                #Temporary
                #if loc <= 9:
                #    continue

                for mut in dict_classifiers[loc]:

                    #Score the mutation
                    prob_ben, prob_neu, prob_del = self.bayes_model_score_combo(dict_classifiers[loc][mut], combo)

                    #Assign a class
                    if prob_ben > prob_neu and prob_ben > prob_del:
                        mut_score = "BEN"
                    elif prob_neu > prob_ben and prob_neu > prob_del:
                        mut_score = "NEU"
                    elif prob_del > prob_ben and prob_del > prob_neu:
                        mut_score = "DEL"

                    #Activity Class
                    class_activity = dict_classifiers[loc][mut][ds_activity]

                    #Stability Class
                    #class_stability = dict_classifiers[loc][mut][ds_stability]

                    #Hardcode this for now:
                    if class_activity == "BEN":
                        list_count_ben.append(mut_score)
                    #elif class_activity == "NEU" and class_stability == "BEN":
                    #    list_count_ben.append(mut_score)      
                    #elif class_activity == "NEU" and class_stability != "BEN":
                    #    list_count_neu.append(mut_score)       
                    elif class_activity == "NEU":
                        list_count_neu.append(mut_score) 
                    elif class_activity == "DEL":
                        list_count_del.append(mut_score)                            
            
            #Tally the mutations
            dict_ben = dict(Counter(list_count_ben))
            dict_neu = dict(Counter(list_count_neu))
            dict_del = dict(Counter(list_count_del))

            #Calculate our rate of del in ben + neu

                
            combo_results = ','.join([

                '0' if 'BEN' not in dict_ben else str(dict_ben['BEN']),
                '0' if 'NEU' not in dict_ben else str(dict_ben['NEU']),
                '0' if 'DEL' not in dict_ben else str(dict_ben['DEL']),
                
                '0' if 'BEN' not in dict_neu else str(dict_neu['BEN']),
                '0' if 'NEU' not in dict_neu else str(dict_neu['NEU']),
                '0' if 'DEL' not in dict_neu else str(dict_neu['DEL']),

                '0' if 'BEN' not in dict_del else str(dict_del['BEN']),
                '0' if 'NEU' not in dict_del else str(dict_del['NEU']),
                '0' if 'DEL' not in dict_del else str(dict_del['DEL']),
                
                ','.join([str(x) for x in combo])
                ]) + '\n'
            print(combo_results)
            output = output + combo_results

        return output


    def mutation_counter(self, list_count_ben, list_count_neu, list_count_del):
        """Count Mutations From BEN, NEU, DEL lists"""

        #Tally the mutations
        dict_ben = dict(Counter(list_count_ben))
        dict_neu = dict(Counter(list_count_neu))
        dict_del = dict(Counter(list_count_del))

        #Return a string with the counts
        return ','.join([

                    '0' if 'BEN' not in dict_ben else str(dict_ben['BEN']),
                    '0' if 'NEU' not in dict_ben else str(dict_ben['NEU']),
                    '0' if 'DEL' not in dict_ben else str(dict_ben['DEL']),
                
                    '0' if 'BEN' not in dict_neu else str(dict_neu['BEN']),
                    '0' if 'NEU' not in dict_neu else str(dict_neu['NEU']),
                    '0' if 'DEL' not in dict_neu else str(dict_neu['DEL']),

                    '0' if 'BEN' not in dict_del else str(dict_del['BEN']),
                    '0' if 'NEU' not in dict_del else str(dict_del['NEU']),
                    '0' if 'DEL' not in dict_del else str(dict_del['DEL']),
                    ])

    def protocol(self):
        """Main entrypoint for the protocol"""
        
        #Create a output log file that we can append to
        with open(self.directory + self.output_prefix + "_" + strftime("%m_%d_%Y") + "-" +
                 strftime("%H_%M_%S") + '_output.txt', 'w') as file_output:
            file_output.write(self.pact_preamble + "\n")

            """
            *****************************************
            Import Classifiers (Required)
            *****************************************
            """
            #Check to see if the section is there
            if not self.obj_cfgparser.has_section('import_classifiers'):           
                print("[Protocols:" + str_protocol_name + " Error] The import_classifiers config file is incorrect.")
                print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [import_classifiers] section.")
                quit()

            #Load the classifiers from the pact file
            dict_classifiers = open_pact_file(self.directory + '/' + self.obj_cfgparser.get('import_classifiers', 'file'))

            """
            *****************************************
            Count Classifiers
            *****************************************
            """
            if self.dict_workflow['bayes_count']:
                self.classifier_count(dict_classifiers)

            """
            *****************************************
            Assess model performance
            *****************************************
            """
            #Define our two datasets
            ds_activity = "amie_ace2_classified"
            ds_activity_fitness = "amie_ace2_fitness"
            #ds_activity = "LGK_Triple_classified"
            #ds_activity_fitness = "LGK_Triple_fitness"
            #ds_stability = "LGK_wt_classified"

            if self.dict_workflow['bayes_model_score']:

                #Setup our counters to measure correct, fp, fn
                list_count_ben = []
                list_count_neu = []
                list_count_del = []

                #Loop each location
                for loc in dict_classifiers:

                    #Temporary
                    #if loc <= 9:
                    #    continue

                    for mut in dict_classifiers[loc]:

                        #Score the mutation
                        prob_ben, prob_neu, prob_del = self.bayes_model_score(dict_classifiers[loc][mut])

                        #Assign a class
                        if prob_ben > prob_neu and prob_ben > prob_del:
                            mut_score = "BEN"
                        elif prob_neu > prob_ben and prob_neu > prob_del:
                            mut_score = "NEU"
                        elif prob_del > prob_ben and prob_del > prob_neu:
                            mut_score = "DEL"

                        #Activity Class
                        class_activity = dict_classifiers[loc][mut][ds_activity]
                        fitness_activity = dict_classifiers[loc][mut][ds_activity_fitness]

                        if fitness_activity == "NaN" or class_activity == "UNCLASSIFIED":
                            continue

                        #Stability Class
                        #class_stability = dict_classifiers[loc][mut][ds_stability]

                        #Hardcode this for now:
                        #if class_activity == "BEN":
                        #    list_count_ben.append(mut_score)
                        #elif class_activity == "NEU" and class_stability == "BEN":
                        #    list_count_ben.append(mut_score)      
                        #elif class_activity == "NEU" and class_stability != "BEN":
                        #    list_count_neu.append(mut_score)       
                        #elif class_activity == "NEU":
                        #    list_count_neu.append(mut_score) 
                        #elif class_activity == "DEL":
                        #    list_count_del.append(mut_score) 

                        #By fitness value
                        if float(fitness_activity) >= -0.3:
                            list_count_ben.append(mut_score)     
                        elif float(fitness_activity) >= -1 and float(fitness_activity) < -0.3:
                            list_count_neu.append(mut_score) 
                        elif float(fitness_activity) < -1:
                            list_count_del.append(mut_score)

                #Print our header
                print("BEN_BEN, BEN_NEU, BEN_DEL, NEU_BEN, NEU_NEU, NEU_DEL, DEL_BEN, DEL_NEU, DEL_DEL")

                #Calculate our rate of del in ben + neu
                print(self.mutation_counter(list_count_ben, list_count_neu, list_count_del))

            if self.dict_workflow['strict_filter_old']:
                #Setup our counters to measure correct, fp, fn
                list_count_ben = []
                list_count_neu = []
                list_count_del = []   

                #Loop each location
                for loc in dict_classifiers:

                    #Temporary
                    #if loc <= 9:
                    #    continue

                    for mut in dict_classifiers[loc]:

                        #Score the mutation
                        mut_score = self.strict_filter_old(dict_classifiers[loc][mut])

                        #Activity Class
                        class_activity = dict_classifiers[loc][mut][ds_activity]
                        fitness_activity = dict_classifiers[loc][mut][ds_activity_fitness]

                        if fitness_activity == "NaN" or class_activity == "UNCLASSIFIED":
                            continue

                        #Stability Class
                        #class_stability = dict_classifiers[loc][mut][ds_stability]

                        #Hardcode this for now:
                        if class_activity == "BEN":
                            list_count_ben.append(mut_score)
                        #elif class_activity == "NEU" and class_stability == "BEN":
                        #    list_count_ben.append(mut_score)      
                        #elif class_activity == "NEU" and class_stability != "BEN":
                        #    list_count_neu.append(mut_score)       
                        elif class_activity == "NEU":
                            list_count_neu.append(mut_score) 
                        elif class_activity == "DEL":
                            list_count_del.append(mut_score)

                        #By fitness value
                        #if float(fitness_activity) >= -0.3:
                        #    list_count_ben.append(mut_score)     
                        #elif float(fitness_activity) >= -1 and float(fitness_activity) < -0.3:
                        #    list_count_neu.append(mut_score) 
                        #elif float(fitness_activity) < -1:
                        #    list_count_del.append(mut_score)
            
                #Print our header
                print("BEN_BEN, BEN_NEU, BEN_DEL, NEU_BEN, NEU_NEU, NEU_DEL, DEL_BEN, DEL_NEU, DEL_DEL")

                #Calculate our rate of del in ben + neu
                print(self.mutation_counter(list_count_ben, list_count_neu, list_count_del))

            if self.dict_workflow['strict_filter_new']:
                #Setup our counters to measure correct, fp, fn
                list_count_ben = []
                list_count_neu = []
                list_count_del = []

                #Loop each location
                for loc in dict_classifiers:

                    #Temporary
                    #if loc <= 9:
                    #    continue

                    for mut in dict_classifiers[loc]:

                        #Score the mutation
                        mut_score = self.strict_filter_new(dict_classifiers[loc][mut])

                        #Activity Class
                        class_activity = dict_classifiers[loc][mut][ds_activity]
                        fitness_activity = dict_classifiers[loc][mut][ds_activity_fitness]

                        if fitness_activity == "NaN" or class_activity == "UNCLASSIFIED":
                            continue

                        #Stability Class
                        #class_stability = dict_classifiers[loc][mut][ds_stability]

                        #Hardcode this for now:
                        if class_activity == "BEN":
                            list_count_ben.append(mut_score)
                        #elif class_activity == "NEU" and class_stability == "BEN":
                        #    list_count_ben.append(mut_score)      
                        #elif class_activity == "NEU" and class_stability != "BEN":
                        #    list_count_neu.append(mut_score)       
                        elif class_activity == "NEU":
                            list_count_neu.append(mut_score) 
                        elif class_activity == "DEL":
                            list_count_del.append(mut_score)
                            
                        #By fitness value
                        #if float(fitness_activity) >= -0.3:
                        #    list_count_ben.append(mut_score)     
                        #elif float(fitness_activity) >= -1 and float(fitness_activity) < -0.3:
                        #    list_count_neu.append(mut_score) 
                        #elif float(fitness_activity) < -1:
                        #    list_count_del.append(mut_score)
            
                #Print our header
                print("BEN_BEN, BEN_NEU, BEN_DEL, NEU_BEN, NEU_NEU, NEU_DEL, DEL_BEN, DEL_NEU, DEL_DEL")

                #Calculate our rate of del in ben + neu
                print(self.mutation_counter(list_count_ben, list_count_neu, list_count_del))

            """
            *****************************************
            Perform a multi-classifier optimization using Bayes probabilities
            *****************************************
            """
            if self.dict_workflow['bayes_combo']:
                file_output.write(self.bayes_combo(dict_classifiers))
        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
