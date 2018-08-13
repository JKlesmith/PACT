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

"""Basal Count - Count the basal rate of classifiers"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Basal Count Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from collections import Counter
from configparser import NoOptionError, NoSectionError
from pact.pact_common import file_checker, pretty_counter_dicts
from pickle import load

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class basal_count:
    """Return the basal classified counts."""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        return

    def basal_count(self, dict_classified):
        """From the config file it loads and combines pact files"""

        #dict_classified is a dict with [location][mutation] = "CLASSIFIED" based on a fitness value
        #Count our classifieds (returns a dict {"CLASSIFED NAME":counts, ...
        dict_counts_screen = dict(Counter(dict_classified[loc][mut]
                    for loc in dict_classified
                    for mut in dict_classified[loc]))

        #Formulate a string of our counts
        str_basal = pretty_counter_dicts(dict_counts_screen)
        print(str_basal)

        return str_basal

if __name__ == '__main__':
    #Remind the user that this file needs to be ran within the context of PACT
    print("[Basal Count Error] This file needs to be ran within the context of PACT.")
