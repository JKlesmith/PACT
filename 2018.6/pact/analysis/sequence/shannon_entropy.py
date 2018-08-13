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

"""shannon_entropy - calculate the shannon entropy for a given residue"""
#Will need to define special codons in the future (if not all are nnk and have 20)

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Shannon Entropy Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

#from ast import literal_eval #for future non nnk codon use
from math import log, pow
from operator import mul
from functools import reduce

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class shannon_entropy_classifier:
    """Calculate the shannon entropy for a given residue with input fitness values."""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Setup the globals"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #General Settings
        try:
            self.directory = self.dict_protocolconfig['directory']
            self.output_prefix = self.config_file.get('global', 'output_prefix')
            self.wtaa = self.config_file.get('global', 'wtaa').upper()
            self.dataset = self.config_file.get('shannon_entropy', 'dataset')
            self.mutation_type = self.config_file.get('shannon_entropy', 'mutation_type').lower()

            if self.mutation_type != "single" or self.mutation_type != "multiple":
                self.mutation_type = "single"

        except NoOptionError:
            print("[Shannon Entropy Error] Missing config file elements.")
            quit()

        return
 
    def import_data(self, dict_merged_datasets):
        """Convert the SSM dict into a flat list at each location"""

        #Empty dict to work with
        dict_enrichments = {}
        
        #Populate our lookup table with the designed location (start 0) and the log2_enrichment value
        for key in dict_merged_datasets[self.dataset]:

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset][key]:

                #Exclude NaN values and stops
                if dict_merged_datasets[self.dataset][key][mut]['log2_enrichment'] == "NaN" or mut == '*':
                    continue

                #Check if the key already exists
                if key in dict_enrichments:
                    dict_enrichments[key] = dict_enrichments[key] + [dict_merged_datasets[self.dataset][key][mut]['log2_enrichment']]
                else:
                    dict_enrichments[key] = [dict_merged_datasets[self.dataset][key][mut]['log2_enrichment']]

        #Currently we have:
        #{238: [-7.46277341387248, -7.015402516210234, -6.975213426094534, -6.457920751911184,
        return dict_enrichments 

    def measured_shannon(self, dict_enrichments):
        """Return the normalized shannon entropy for a design region"""

        #Empty dict to work with
        dict_col_sum = {}
        dict_col_counts = {}
        dict_col_pvals = {}
        dict_col_se = {}
        dict_col_se_total = {}
        dict_total_mutations = {}
        dict_normalized_se = {}

        #Populate our lookup table with the designed location (start 1) and the log2
        for key in dict_enrichments:
           
            #Get the location column sum (2^enrichment)
            dict_col_sum[key] = sum(list(pow(2, e) for e in dict_enrichments[key]))
            
            #Get the number of values in each column (+1 for wild-type)
            dict_col_counts[key] = len(list(pow(2, e) for e in dict_enrichments[key])) + 1

            #Get the pvals for each mutation (2^enrichment/column sum)
            dict_col_pvals[key] = list((pow(2, e)/dict_col_sum[key]) 
                                       for e in dict_enrichments[key])

            #Calculate the entropy value per mutation -1 * pval * log2(pval) (WT is already included from above)
            dict_col_se[key] = list(-1 * e * log(e, 2) for e in dict_col_pvals[key])

            #Get the column total entropy
            dict_col_se_total[key] = sum(e for e in dict_col_se[key])

            #Get the total possible mutations per site (20*20...)
            if self.mutation_type == "single":
                dict_total_mutations[key] = 20
            else:
                dict_total_mutations[key] = reduce(mul, list(e for e in len(key.split(',')) * [20]))

            #Return the normalized column total se (se total * (log2(total muts possible)/log2(muts actually seen))
            dict_normalized_se[key] = dict_col_se_total[key] * (log(dict_total_mutations[key], 2)/log(dict_col_counts[key], 2))

        return dict_normalized_se
    
    def theoretical_shannon(self, dict_enrichments):
        """Return the theoretical shannon entropy for a given design site."""

        #Empty dict to work with
        dict_col_sum = {}
        dict_col_pvals = {}
        dict_col_se = {}
        dict_col_se_total = {}
        dict_total_mutations = {}

        #Populate our lookup table with the designed location (start 1) and the log2
        for key in dict_enrichments:

            #Get the total possible mutations per site (20*20...) (1*20 for single libraries)
            if self.mutation_type == "single":
                dict_total_mutations[key] = 20
            else:
                dict_total_mutations[key] = reduce(mul, list(e for e in len(key.split(',')) * [20]))

            #Get the location column sum (2^enrichment)
            dict_col_sum[key] = sum(list(pow(2, e) for e in dict_total_mutations[key] * [1]))

            #Get the pvals for each mutation (2^enrichment/column sum)
            dict_col_pvals[key] = list((pow(2, e)/dict_col_sum[key]) for e in dict_total_mutations[key] * [1])

            #Calculate the entropy value per mutation -1 * pval * log2(pval)
            dict_col_se[key] = list(-1 * e * log(e, 2) for e in dict_col_pvals[key])

            #Get the column total entropy
            dict_col_se_total[key] = sum(e for e in dict_col_se[key])
       
        return dict_col_se_total

    def shannon_entropy(self, dict_merged_datasets):
        """Main entrypoint for the class"""

        #Convert our fitness dataset into a list at each site
        dict_fitness_lists = self.import_data(dict_merged_datasets)

        #Calculate the measured shannon entropy
        dict_measured_se = self.measured_shannon(dict_fitness_lists)

        #Calculate the theoretical shannon entropy
        dict_theo_se = self.theoretical_shannon(dict_fitness_lists)

        #Compare the two dicts
        dict_results = {
            key : [round(dict_measured_se[key], 4), round(dict_theo_se[key], 4), round((dict_measured_se[key]/dict_theo_se[key]), 4)]
            for key in dict_measured_se}

        #Output our data to the screen
        print("Location", "Measured SE", "Theoretical SE", "Ratio Measured/Theoretical")
        print('\n'.join([
            str(key) + ',' + str(dict_results[key][0])
            + ',' + str(dict_results[key][1])
            + ',' + str(dict_results[key][2])
            for key in dict_results
            ]))
        
        #Open the output file and write it
        with open(self.directory + self.output_prefix + "_shannon.csv", "w") as csv_out:
            csv_out.write("Location,Measured SE,Theoretical SE,Ratio Measured/Theoretical\n")
            csv_out.write('\n'.join([
                str(key) + ',' + str(dict_results[key][0])
                + ',' + str(dict_results[key][1])
                + ',' + str(dict_results[key][2])
                for key in dict_results
                ]))
        
        return dict_results

if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[Shannon Entropy Error] This classifier needs to be ran within the context of PACT.")
