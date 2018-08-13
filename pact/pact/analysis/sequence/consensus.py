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

"""consensus - functions to probe back to consensus mutations"""

from sys import version_info

#Setup our classifier name
str_classifier_name = "Consensus"

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[" + str_classifier_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError, NoSectionError
from pact_common import create_csv, get_bool

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class consensus:
    """Functions to probe back to consensus mutations"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""

        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Classifier specific variables
        try:
            self.wtaa = self.config_file.get('global', 'wtaa')
            self.directory = self.dict_protocolconfig['directory']
            self.output_prefix = self.config_file.get('global', 'output_prefix')
        except NoOptionError:
            print("[" + str_classifier_name + " Error] Missing config file elements.")
            quit()

        return

    def site_max(self, dict_pssm, loc, column):
        """Get the site max value from the pssm dict"""
        return max([int(dict_pssm[loc][mut][column]) for mut in dict_pssm[loc]])  
    
    def site_count(self, dict_pssm, loc, column):
        """Get the count of conserved values per site"""

        #Column zero is PSSM
        if column == 0:
            return len([1 for mut in dict_pssm[loc] if int(dict_pssm[loc][mut][column]) >= 0])
        else:
            return len([1 for mut in dict_pssm[loc] if int(dict_pssm[loc][mut][column]) > 0])

    def wt_consensus(self, dict_pssm):
        """Return information on the wild-type sequence"""

        #Print the base wt stats
        str_header = ",".join([
            "Location",
            "WT Residue",
            "PSSM Value",
            "Percent Observed",
            "Max Site PSSM",
            "Max Site Percent",
            "Site PSSM Above Zero",
            "Site Percent Non-Zero",
            "PSSM WT Max Cons?",
            "Percent WT Max Cons?",
            ])

        str_body = ""

        #Create a dict to return the info
        dict_wtcons = {}

        #At each site return the % observed for the wild-type sequence
        for i in range(0, len(self.wtaa)):

            #Make the site
            dict_wtcons[i + 1] = {}

            #Store the wt pssm and percent value
            wt_pssm = int(dict_pssm[i + 1][self.wtaa[i]][0])
            wt_percent = int(dict_pssm[i + 1][self.wtaa[i]][1])

            #Store the max and total conserved site-wise counts
            max_pssm = self.site_max(dict_pssm, i + 1, 0)
            max_percent = self.site_max(dict_pssm, i + 1, 1)
            pssm_cons_count = self.site_count(dict_pssm, i + 1, 0)
            perc_cons_count = self.site_count(dict_pssm, i + 1, 1)

            #Is the wt sequence the site max?
            if wt_pssm == max_pssm:
                wt_max_pssm = 1
            else:
                wt_max_pssm = 0

            if wt_percent == max_percent:
                wt_max_percent = 1
            else:
                wt_max_percent = 0

            str_body = str_body + ",".join(map(str, [
                i + 1,
                self.wtaa[i],
                wt_pssm,
                wt_percent,
                max_pssm,
                max_percent,
                pssm_cons_count,
                perc_cons_count,
                wt_max_pssm,
                wt_max_percent,
                "\n"
                ]))

            dict_wtcons[i + 1] = {
                'wt_resi':self.wtaa[i],
                'wt_pssm':wt_pssm,
                'wt_percent':wt_percent,
                'max_pssm':max_pssm,
                'max_percent':max_percent,
                'pssm_cons_count':pssm_cons_count,
                'percent_cons_count':perc_cons_count,
                'wt_max_pssm':wt_max_pssm,
                'wt_max_percent':wt_max_percent
                }

        #Output our information
        print(str_header)
        print(str_body)

        #Write CSV
        print(create_csv(self.directory + self.output_prefix + "_wt_consensus.csv", str_header, str_body))

        return dict_wtcons

    def wtcons_count_class(self, dict_wtcons, dict_classified, dataset):
        """What is the probability of finding a classifed mutation per wt consensus"""

        #Get the unique class types by converting to a set
        set_class_types = set([dict_classified[loc][mut] 
                               for loc in dict_classified 
                               for mut in dict_classified[loc]])
        
        #Loop each class type
        for class_type in set_class_types:

            #Setup our list
            list_wtcons_pssm = []
            list_wtcons_percent = []

            #Loop the locations
            for loc in dict_classified:

                #Loop the mutations
                for mut in dict_classified[loc]:

                    #Only record from our class type
                    if dict_classified[loc][mut] != class_type:
                        continue

                    #Skip stops
                    if mut == "*":
                        continue

                    #Skip WT
                    if mut == dict_wtcons[loc]['wt_resi']:
                        continue

                    #Append our WT data
                    list_wtcons_pssm.append(dict_wtcons[loc]['wt_pssm'])
                    list_wtcons_percent.append(dict_wtcons[loc]['wt_percent'])

            #Create our output strings
            str_header = "WT PSSM,WT Percent"
            str_body = '\n'.join([str(x[0]) + "," + str(x[1]) for x in zip(list_wtcons_pssm, list_wtcons_percent)])

            print("Class: " + class_type)
            print(str_header)
            print(str_body)

            #Write CSV
            print(create_csv(self.directory + self.output_prefix + "_Dataset_" + dataset + "_Class_" +
                            class_type + "_wt_consensus.csv", str_header, str_body))

        return

    def nonconserved_sites(self, dict_wtcons, dict_pssm, dict_classified, dataset):
        """What is the classification at sites where WT is not conserved if we mutate to the top/any conserved"""

        #Get the unique class types by converting to a set
        set_class_types = set([dict_classified[loc][mut] 
                               for loc in dict_classified 
                               for mut in dict_classified[loc]])

        str_output = ""
        
        #Loop each class type
        for class_type in set_class_types:

            #Setup our list
            list_top_pssm = []
            list_top_percent = []
            list_any_pssm = []
            list_any_percent = []
            list_basal_pssm = []
            list_basal_percent = []

            #Loop the locations
            for loc in dict_classified:

                #Loop the mutations
                for mut in dict_classified[loc]:

                    #Only record from our class type
                    if dict_classified[loc][mut] != class_type:
                        continue

                    #Skip stops
                    if mut == "*":
                        continue

                    #Skip WT
                    if mut == dict_wtcons[loc]['wt_resi']:
                        continue

                    #Append our classifier data
                    if (dict_wtcons[loc]['wt_max_pssm'] == 0 and 
                        int(dict_pssm[loc][mut][0]) == dict_wtcons[loc]['max_pssm']):
                        list_top_pssm.append(1)

                    if (dict_wtcons[loc]['wt_max_percent'] == 0 and
                        int(dict_pssm[loc][mut][1]) == dict_wtcons[loc]['max_percent']):
                        list_top_percent.append(1)

                    if (dict_wtcons[loc]['wt_max_pssm'] == 0 and 
                        int(dict_pssm[loc][mut][0]) >= 0):
                        list_any_pssm.append(1)

                    if (dict_wtcons[loc]['wt_max_percent'] == 0 and
                        int(dict_pssm[loc][mut][1]) > 0):
                        list_any_percent.append(1)

                    if int(dict_pssm[loc][mut][0]) >= 0:
                        list_basal_pssm.append(1)

                    if int(dict_pssm[loc][mut][1]) > 0:
                        list_basal_percent.append(1)
            
            str_output = str_output + ",".join([
                dataset,
                class_type,
                str(len(list_top_pssm)),
                str(len(list_top_percent)),                
                str(len(list_any_pssm)),
                str(len(list_any_percent)),
                str(len(list_basal_pssm)),
                str(len(list_basal_percent)),
                ]) + "\n"
 
        #Setup the output header
        str_header = ",".join([
            "Dataset",
            "Classification",
            "Top PSSM",
            "Top Percent",              
            "Any PSSM",
            "Any Percent",
            "Basal PSSM",
            "Basal Percent"
            ])

        print(str_header)
        print(str_output)

        return str_header + "\n" + str_output

    def nonconserved_mutations(self, dict_wtcons, dict_pssm, dict_classified, dataset):
        """What is the classification at sites where WT is not conserved if we mutate to non conserved residue"""

        #Get the unique class types by converting to a set
        set_class_types = set([dict_classified[loc][mut] 
                               for loc in dict_classified 
                               for mut in dict_classified[loc]])

        str_output = ""
        
        #Loop each class type
        for class_type in set_class_types:

            #Setup our list
            list_wtnc_pssm = []
            list_wtnc_percent = []
            list_basal_pssm = []
            list_basal_percent = []

            #Loop the locations
            for loc in dict_classified:

                #Loop the mutations
                for mut in dict_classified[loc]:

                    #Only record from our class type
                    if dict_classified[loc][mut] != class_type:
                        continue

                    #Skip stops
                    if mut == "*":
                        continue

                    #Skip WT
                    if mut == dict_wtcons[loc]['wt_resi']:
                        continue

                    #Append our classifier data
                    if (dict_wtcons[loc]['wt_max_pssm'] == 0 and 
                        int(dict_pssm[loc][mut][0]) < 0):
                        list_wtnc_pssm.append(1)

                    if (dict_wtcons[loc]['wt_max_percent'] == 0 and
                        int(dict_pssm[loc][mut][1]) == 0):
                        list_wtnc_percent.append(1)

                    if int(dict_pssm[loc][mut][0]) < 0:
                        list_basal_pssm.append(1)

                    if int(dict_pssm[loc][mut][1]) == 0:
                        list_basal_percent.append(1)
            
            str_output = str_output + ",".join([
                dataset,
                class_type,
                str(len(list_wtnc_pssm)),
                str(len(list_wtnc_percent)),                
                str(len(list_basal_pssm)),
                str(len(list_basal_percent)),
                ]) + "\n"
 
        #Setup the output header
        str_header = ",".join([
            "Dataset",
            "Classification",             
            "WTNC PSSM",
            "WTNC Percent",
            "Basal PSSM",
            "Basal Percent"
            ])

        print(str_header)
        print(str_output)

        return str_header + "\n" + str_output

    def cons_count_setvset(self, dict_wtcons, dict_pssm, dict_classified):
        """Report the set vs set class v class consensus for the wt residue and the mutation"""

        #Classifier specific variables
        try:
            self.datasetx = self.config_file.get("consensus", "dataset_x")
            self.datasety = self.config_file.get("consensus", "dataset_y")
        except NoOptionError:
            print("[" + str_classifier_name + " Error] Missing config file elements.")
            quit()

        #Get the unique class types by converting to a set
        set_class_types = set([dict_classified[dataset][loc][mut] 
                               for dataset in dict_classified
                               for loc in dict_classified[dataset] 
                               for mut in dict_classified[dataset][loc]])
        
        #Loop each class type
        for class_type_x in set_class_types:
            for class_type_y in set_class_types:

                #Setup our list
                list_wtcons_pssm = []
                list_wtcons_percent = []
                list_mutcons_pssm = []
                list_mutcons_percent = []

                #Loop the locations
                for loc in dict_classified[self.datasetx]:

                    #Loop the mutations
                    for mut in dict_classified[self.datasetx][loc]:

                        #Only record from our class type
                        if dict_classified[self.datasetx][loc][mut] != class_type_x:
                            continue

                        #Only record from our class type
                        if dict_classified[self.datasety][loc][mut] != class_type_y:
                            continue

                        #Skip stops
                        if mut == "*":
                            continue

                        #Skip WT
                        if mut == dict_wtcons[loc]['wt_resi']:
                            continue

                        #Append our WT data
                        list_wtcons_pssm.append(dict_wtcons[loc]['wt_pssm'])
                        list_wtcons_percent.append(dict_wtcons[loc]['wt_percent'])
                        list_mutcons_pssm.append(int(dict_pssm[loc][mut][0]))
                        list_mutcons_percent.append(int(dict_pssm[loc][mut][1]))

                #Create our output strings
                str_header = "WT PSSM,WT Percent,Mut PSSM,Mut Percent"
                str_body = '\n'.join([str(x[0]) + "," + str(x[1]) + "," + str(x[2]) + "," + str(x[3]) 
                                      for x in zip(list_wtcons_pssm, list_wtcons_percent, list_mutcons_pssm, list_mutcons_percent)])
            
                print("Dataset X: " + self.datasetx + "Dataset Y: " + self.datasety)
                print("Class X: " + class_type_x + "Class Y: " + class_type_y)
                print(str_header)
                print(str_body)

                #Write CSV
                print(create_csv(self.directory + self.output_prefix + "_" + self.datasetx + "_" +
                                class_type_x + "_" + self.datasety + "_" + class_type_y + 
                                "_setvset_cons.csv", str_header, str_body))

        return

    def nonconserved_sites_burial(self, dict_wtcons, dict_pssm, dict_classified, dataset, dict_pdb, equality):
        """What is the classification at sites where WT is not conserved if we mutate to the top/any conserved"""

        #Get the unique class types by converting to a set
        set_class_types = set([dict_classified[loc][mut] 
                               for loc in dict_classified 
                               for mut in dict_classified[loc]])

        str_output = ""

        try:
            pdb_file = self.config_file.get("consensus", "pdb_file")
            frac_burial = float(self.config_file.get("consensus", "frac_burial"))
            chain = self.config_file.get("consensus", "chain").upper()
        except NoOptionError:
            print("[" + str_classifier_name + " Error] Missing config file elements.")
            quit()
        
        #Loop each class type
        for class_type in set_class_types:

            #Setup our list
            list_top_pssm = []
            list_top_percent = []
            list_any_pssm = []
            list_any_percent = []
            list_basal_pssm = []
            list_basal_percent = []

            #Loop the locations
            for loc in dict_classified:

                #Loop the mutations
                for mut in dict_classified[loc]:

                    #Only record from our class type
                    if dict_classified[loc][mut] != class_type:
                        continue

                    #Skip stops
                    if mut == "*":
                        continue

                    #Skip WT
                    if mut == dict_wtcons[loc]['wt_resi']:
                        continue

                    #Check for burial
                    if get_bool(dict_pdb[pdb_file]['dssp'][chain][loc]['frac_burial'], equality, frac_burial):
                        continue

                    #Append our classifier data
                    if (dict_wtcons[loc]['wt_max_pssm'] == 0 and 
                        int(dict_pssm[loc][mut][0]) == dict_wtcons[loc]['max_pssm']):
                        list_top_pssm.append(1)

                    if (dict_wtcons[loc]['wt_max_percent'] == 0 and
                        int(dict_pssm[loc][mut][1]) == dict_wtcons[loc]['max_percent']):
                        list_top_percent.append(1)

                    if (dict_wtcons[loc]['wt_max_pssm'] == 0 and 
                        int(dict_pssm[loc][mut][0]) >= 0):
                        list_any_pssm.append(1)

                    if (dict_wtcons[loc]['wt_max_percent'] == 0 and
                        int(dict_pssm[loc][mut][1]) > 0):
                        list_any_percent.append(1)

                    if int(dict_pssm[loc][mut][0]) >= 0:
                        list_basal_pssm.append(1)

                    if int(dict_pssm[loc][mut][1]) > 0:
                        list_basal_percent.append(1)
            
            str_output = str_output + ",".join([
                dataset,
                class_type,
                str(len(list_top_pssm)),
                str(len(list_top_percent)),                
                str(len(list_any_pssm)),
                str(len(list_any_percent)),
                str(len(list_basal_pssm)),
                str(len(list_basal_percent)),
                ]) + "\n"
 
        #Setup the output header
        str_header = ",".join([
            "Dataset",
            "Classification",
            "Top PSSM",
            "Top Percent",              
            "Any PSSM",
            "Any Percent",
            "Basal PSSM",
            "Basal Percent"
            ])

        print(str_header)
        print(str_output)

        return str_header + "\n" + str_output

if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[" + str_classifier_name + " Error] This classifier needs to be ran within the context of PACT.") 
