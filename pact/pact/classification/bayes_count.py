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
    print("[Bayes Count Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from collections import Counter
from pact.pact_common import get_bool

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class bayes_count:
    """Return the basal classified counts."""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Get the config file items
        self.class_key = self.config_file.get('bayes_count', 'classification_key')
        print("Classified dataset: " + self.class_key)

        return

    def classifier_count_list(self, list_ben, list_neu, list_del, dict_key1, list_group1):
        """Count classifiers that are numberic vs group of items in list (ranges vs list of items)"""

        #List of dicts
        list_ben_count = 0
        list_neu_count = 0
        list_del_count = 0

        #BEN List
        for dict_mutation in list_ben:
            if dict_mutation[dict_key1] in list_group1:
                list_ben_count += 1

        #NEU List
        for dict_mutation in list_neu:
            if dict_mutation[dict_key1] in list_group1:
                list_neu_count += 1 

        #NEU List
        for dict_mutation in list_del:
            if dict_mutation[dict_key1] in list_group1:
                list_del_count += 1

        #Print
        print("BEN: " + str(list_ben_count) + " NEU: " + str(list_neu_count) + " DEL: " + str(list_del_count))
        return

    def classifier_count_listlist(self, list_ben, list_neu, list_del, dict_key1, dict_key2, 
                                   list_group1, list_group2):
        """Count classifiers that are numberic vs group of items in list (ranges vs list of items)"""

        #List of dicts
        list_ben_count = 0
        list_neu_count = 0
        list_del_count = 0

        #BEN List
        for dict_mutation in list_ben:
            if dict_mutation[dict_key1] in list_group1 and dict_mutation[dict_key2] in list_group2:
                list_ben_count += 1

        #NEU List
        for dict_mutation in list_neu:
            if dict_mutation[dict_key1] in list_group1 and dict_mutation[dict_key2] in list_group2:
                list_neu_count += 1 

        #NEU List
        for dict_mutation in list_del:
            if dict_mutation[dict_key1] in list_group1 and dict_mutation[dict_key2] in list_group2:
                list_del_count += 1

        #Print
        print("BEN: " + str(list_ben_count) + " NEU: " + str(list_neu_count) + " DEL: " + str(list_del_count))
        return

    def classifier_count_rangelist(self, list_ben, list_neu, list_del, dict_key1, dict_key2, 
                                   range_start, relation_start, range_end, relation_end, list_group2):
        """Count classifiers that are numberic vs group of items in list (ranges vs list of items)"""

        #List of dicts
        list_ben_count = 0
        list_neu_count = 0
        list_del_count = 0

        #BEN List
        for dict_mutation in list_ben:
            if (get_bool(float(dict_mutation[dict_key1]), relation_start, range_start) and 
                get_bool(float(dict_mutation[dict_key1]), relation_end, range_end) and
                dict_mutation[dict_key2] in list_group2):
                list_ben_count += 1

        #NEU List
        for dict_mutation in list_neu:
            if (get_bool(float(dict_mutation[dict_key1]), relation_start, range_start) and 
                get_bool(float(dict_mutation[dict_key1]), relation_end, range_end) and
                dict_mutation[dict_key2] in list_group2):
                list_neu_count += 1 

        #NEU List
        for dict_mutation in list_del:
            if (get_bool(float(dict_mutation[dict_key1]), relation_start, range_start) and 
                get_bool(float(dict_mutation[dict_key1]), relation_end, range_end) and
                dict_mutation[dict_key2] in list_group2):
                list_del_count += 1

        #Print
        print("BEN: " + str(list_ben_count) + " NEU: " + str(list_neu_count) + " DEL: " + str(list_del_count))
        return

    def classifier_count_keycompare(self, list_ben, list_neu, list_del, dict_key1, dict_key2, relation):
        """Count classifiers that are numberic (ranges vs ranges)"""

        #List of dicts
        list_ben_count = 0
        list_neu_count = 0
        list_del_count = 0

        #BEN List
        for dict_mutation in list_ben:
            if get_bool(float(dict_mutation[dict_key1]), relation, float(dict_mutation[dict_key2])):
                list_ben_count += 1

        #NEU List
        for dict_mutation in list_neu:
            if get_bool(float(dict_mutation[dict_key1]), relation, float(dict_mutation[dict_key2])):
                list_neu_count += 1 

        #NEU List
        for dict_mutation in list_del:
            if get_bool(float(dict_mutation[dict_key1]), relation, float(dict_mutation[dict_key2])):
                list_del_count += 1

        #Print
        print("BEN: " + str(list_ben_count) + " NEU: " + str(list_neu_count) + " DEL: " + str(list_del_count))
        return

    def classifier_count_range(self, list_ben, list_neu, list_del, dict_key, 
                               range_start, relation_start, range_end, relation_end):
        """Count classifiers that are numberic (ranges)"""

        #List of dicts
        list_ben_count = 0
        list_neu_count = 0
        list_del_count = 0

        #BEN List
        for dict_mutation in list_ben:
            if (get_bool(float(dict_mutation[dict_key]), relation_start, range_start) and 
                get_bool(float(dict_mutation[dict_key]), relation_end, range_end)):
                list_ben_count += 1

        #NEU List
        for dict_mutation in list_neu:
            if (get_bool(float(dict_mutation[dict_key]), relation_start, range_start) and 
                get_bool(float(dict_mutation[dict_key]), relation_end, range_end)):
                list_neu_count += 1 

        #NEU List
        for dict_mutation in list_del:
            if (get_bool(float(dict_mutation[dict_key]), relation_start, range_start) and 
                get_bool(float(dict_mutation[dict_key]), relation_end, range_end)):
                list_del_count += 1

        #Print
        print("BEN: " + str(list_ben_count) + " NEU: " + str(list_neu_count) + " DEL: " + str(list_del_count))
        return

    def bayes_count(self, dict_classifiers):
        """From the config file it loads and combines pact files"""

        # Temp variables and settings
        # BEN: LGK_Triple_classified:BEN or (LGK_wt_classified:BEN and LGK_Triple_classified:NEU)
        # NEU: LGK_Triple_classified:NEU where LGK_wt_classified:NEU or LGK_wt_classified:DEL
        # DEL: LGK_Triple_classified:DEL

        #Define our two datasets
        ds_activity = "LGK_Triple_classified"
        ds_stability = "LGK_wt_classified"
        
        #List of dicts
        list_ben = []
        list_neu = []
        list_del = []

        #Split our main dict into the three lists so it will be easier to parse it
        for loc in dict_classifiers:

            #Skip 9 and below
            if loc <= 9:
                continue

            #Loop the mutations
            for mut in dict_classifiers[loc]:

                #Activity Class
                class_activity = dict_classifiers[loc][mut][ds_activity]

                #Stability Class
                class_stability = dict_classifiers[loc][mut][ds_stability]

                #Hardcode this for now:
                if class_activity == "BEN":
                    list_ben.append(dict_classifiers[loc][mut])
                elif class_activity == "NEU" and class_stability == "BEN":
                    list_ben.append(dict_classifiers[loc][mut])
                elif class_activity == "NEU" and class_stability != "BEN":
                    list_neu.append(dict_classifiers[loc][mut])
                elif class_activity == "DEL":
                    list_del.append(dict_classifiers[loc][mut])
        
        #Count the number of items
        print("BEN: " + str(len(list_ben)) + " NEU: " + str(len(list_neu)) + " DEL: " + str(len(list_del)))

        return list_ben, list_neu, list_del

if __name__ == '__main__':
    #Remind the user that this file needs to be ran within the context of PACT
    print("[Bayes Count Error] This file needs to be ran within the context of PACT.")
