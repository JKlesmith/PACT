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

"""Calculate the fitness of mutations"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Fitness Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from configparser import NoOptionError
from pickle import load
from math import log, sqrt, exp, pow, e, pi, factorial
from statistics import stdev, median, mean, StatisticsError
from pact.pact_common import file_checker, save_pact_file

#Setup for FACS fitness metric
try:
    from scipy.special import erfinv
    erfinv_import = True
except ImportError:
    print("[Fitness Error] SciPy not installed, trying alt. implementation of erfinv.")
    try:
        from pact.external_scripts.pyerf import erfinv
        erfinv_import = True
    except ImportError:
        print("[Fitness Error] Import of pyerf 1.0.1 for erfinv failed.")
        erfinv_import = False

#Setup for the ttest
try:
    from scipy.stats import t
    sf_import = True
except ImportError:
    print("[Fitness Error] SciPy not installed, p-values cannot be calculated.")
    sf_import = False

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class fitness:
    """Calculate the fitness of mutations"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("fitness"):
            print("[Fitness Error] the config file is missing the section [fitness]")
            quit()

        #Set the class specific variables assigned by the protocol
        self.directory = self.dict_protocolconfig['directory']
        self.out_prefix = self.dict_protocolconfig['Out_Prefix']
        self.wtdna = self.dict_protocolconfig['WTDNA']
        self.wtaa = self.dict_protocolconfig['WTAA']
        self.FirstAAMutated = int(self.dict_protocolconfig['FirstAAMutated'])-1
        self.LastAAMutated = int(self.dict_protocolconfig['LastAAMutated'])-1
        self.WTAARegion = self.wtaa[self.FirstAAMutated:self.LastAAMutated + 1]
        self.WTDNARegion = self.wtdna[self.FirstAAMutated*3:(self.LastAAMutated + 1) * 3]
        
        #Import the library mode (Single or Multiple)
        if self.dict_protocolconfig['library_type'].lower() == 'single':
            self.library_type = "single"
        elif self.dict_protocolconfig['library_type'].lower() == 'multiple':
            self.library_type = "multiple"
        else:
            print("[Fitness Error] Unknown mode (expect single or multiple).")
            quit()

        #Import settings from the config file on file imports
        try:
            #Test if we have a special input file
            if len(self.config_file.get('fitness', 'pact_enrichment_summary')) > 0:
                file_pact_summary = self.config_file.get('fitness', 'pact_enrichment_summary')
            else:
                file_pact_summary = self.out_prefix + "_enrichment_summary.pact"

            #Test if we have a special input file
            if len(self.config_file.get('fitness', 'pact_enrichment_accept_nonsynon')) > 0:
                file_pact_accepted = self.config_file.get('fitness', 'pact_enrichment_accept_nonsynon')
            else:
                file_pact_accepted = self.out_prefix + "_enrichment_accept_nonsynon.pact"

            #Test if we have a special input file
            if len(self.config_file.get('fitness', 'pact_enrichment_wtsynon')) > 0:
                file_pact_wtsynon = self.config_file.get('fitness', 'pact_enrichment_wtsynon')
            else:
                file_pact_wtsynon = self.out_prefix + "_enrichment_wtsynon.pact"

            #Check our inported files and load
            if file_checker(file_pact_summary):
                with open(file_pact_summary, 'rb') as file_summary:
                    self.dict_summary = load(file_summary)

            if file_checker(file_pact_accepted):
                with open(file_pact_accepted, 'rb') as file_accepted:
                    self.dict_accepted = load(file_accepted)        

            if file_checker(file_pact_wtsynon):
                with open(file_pact_wtsynon, 'rb') as file_wtsynon:
                    self.dict_wtsynon = load(file_wtsynon)

        except NoOptionError:
            print("[Fitness Error] The fitness config file is incorrect.")
            print("[Fitness Error] Missing an option flag that starts with file_.")
            quit()

        #Check to see if we have a manual log2 enrichment
        if len(self.config_file.get('fitness', 'manual_log2')):
            try:
                self.wt_log2 = float(self.config_file.get('fitness', 'manual_log2'))
            except:
                self.wt_log2 = self.dict_summary['log2_wildtype']
        else:
            self.wt_log2 = self.dict_summary['log2_wildtype']

        #Handle the import of the fitness metric
        try:
            #Import the fitness metric type
            self.metric = self.config_file.get('fitness', 'metric').lower()

            #Check to see if the metric is any that we accept
            if self.metric not in ['e-wt', 'growth', 'facs']:
                print("[Fitness Error] Unknown fitness metric type.")
                quit()

            #Check to see if we have a number of doublings
            if self.metric == "growth":
                try:
                    self.growth_gp = float(self.config_file.get('fitness', 'growth_gp'))
                    self.growth_wildtype = ((self.wt_log2 / self.growth_gp) + 1)
                except:
                    print("[Fitness Error] The fitness config file is incorrect for growth_gp.")
                    quit()

            #Check to see if we have our sd or %collected for facs
            if self.metric == "facs":
                #See if erfinv is imported
                if not erfinv_import:
                    print("[Fitness Error] Inverse error function can not be loaded for FACS.")
                    quit()
            
                try:
                    self.facs_sd = float(self.config_file.get('fitness', 'facs_sd')) #Standard Deviation
                    self.facs_pc = float(self.config_file.get('fitness', 'facs_pc')) #Percent collected

                    #Print info from config file
                    print("[Fitness] FACS SD: " + str(round(self.facs_sd, 3)))
                    print("[Fitness] FACS percent collected: " + str(round(self.facs_pc, 3)))
                    print("[Fitness] WT log2 enrichment: " + str(round(self.wt_log2, 3)))

                    #Calculate the inner value
                    wt_facs_inner = 1 - self.facs_pc * pow(2, (self.wt_log2 + 1))

                    #The erfinv is -1 to 1 inclusive
                    if wt_facs_inner <= -1:
                        print("[Fitness] the WT facs equation is less than or equal to -1 for erfinv")
                        print("[Fitness] setting to -0.9999")
                        wt_facs_inner = -0.9999
                    elif wt_facs_inner >= 1:
                        print("[Fitness] the WT facs equation is greater than or equal to 1 for erfinv")
                        print("[Fitness] setting to 0.9999")
                        wt_facs_inner = 0.9999
            
                    #Theoretical maximum enrichment
                    self.e_max_theo_facs = -log(self.facs_pc, 2)
                    self.e_wildtype_facs = erfinv(wt_facs_inner)

                    #Print Theoretical Values
                    print("[Fitness] FACS equation max theo. enrichment: " + str(round(self.e_max_theo_facs, 3)))
                    print("[Fitness] FACS equation wild-type: " + str(round(self.e_wildtype_facs, 3)))

                except TypeError:
                    print("[Fitness Error] The fitness config file is incorrect for facs_sd or facs_pc.")
                    quit()

        except NoOptionError:
            print("[Fitness Error] The fitness config file is incorrect.")
            print("[Fitness Error] There is something wrong with the name of a option flag.")
            quit()

        #Check to see what type of expect value to calculate
        try:
            self.evalue_type = self.config_file.get('fitness', 'evalue_type').lower()

            #Check to see if the metric is any that we accept
            if self.evalue_type not in ['growth', 'facs']:
                self.evalue_type = None
                print("[Fitness] Not performing a t-test for variant significance.")

            #Prepare growth
            if self.evalue_type == "growth":
                #Create the pascal triangle list
                self.list_pascal = [[int((factorial(row)) / ((factorial(k)) * factorial(row - k))) 
                                 for k in range(row + 1)] 
                                 for row in range(int(self.growth_gp) + 1)]

            #Prepare FACS
            if self.evalue_type == "facs":
                try:
                    self.evalue_facs_cellcount = int(self.config_file.get('fitness', 'evalue_facs_cellcount'))
                except:
                    self.evalue_type = None
        except:
            self.evalue_type = None
            print("[Fitness] Not performing a t-test for variant significance.")

        #Import the mutation design list
        self.list_mutation_design = []
        for group in literal_eval(self.dict_protocolconfig['mutcodons']):
            #If we have n between two numbers
            if 'n' in group:
                #We need exactly three to expand the range
                if len(group) == 3:
                    #Get the start and end position, then iterate through the two points
                    self.list_mutation_design.append(list(i for i in range(group[0], group[-1] + 1)))
                else:
                    print("[Fitness Error] The codon list is missing either a start or end point around the 'n' location.")
                    quit()
            else:
                #Add the other ranges without changing them
                self.list_mutation_design.append(group)
        return

    def fitness_metric_growth(self, ei_mutation):
        """Return the metric for growth"""

        mutant = ((ei_mutation / self.growth_gp) + 1)
        
        #Return NaN if we have a very negative value
        if (mutant / self.growth_wildtype) < 0:
            return "NaN"
        else:
            return log(mutant / self.growth_wildtype, 2)

    def fitness_metric_facs(self, ei_mutation):
        """Return the metric for FACS"""

        #Handle the theoretical max enrichment
        try:
            if ei_mutation >= self.e_max_theo_facs:
                mutant = erfinv(1 - self.facs_pc * pow(2, ((self.e_max_theo_facs - 0.001) + 1)))    
            else:
                mutant = erfinv(1 - self.facs_pc * pow(2, (ei_mutation + 1)))
        except:
            print("[Fitness Error] Error within the FACS normalization equation.")
            print("[Fitness Error] ei_mutation: " + str(ei_mutation))
            print("[Fitness Error] theoretical max enrichment: " + str(self.e_max_theo_facs))
            print("[Fitness Error] wild-type enrichment: " + str(self.e_wildtype_facs))
            print("[Fitness Error] percent collected: " + str(self.facs_pc))
            print("[Fitness Error] SD: " + str(self.facs_sd))
            print("[Fitness Error] mutant calculated: " + str(mutant))
        
        return log(e, 2) * sqrt(2) * self.facs_sd * (self.e_wildtype_facs - mutant)

    def welch_pval(self, n1, n2, s1, s2, x1, x2):
        """Calculate the p-value for a Welch ttest"""

        #Calculate the t value
        t_value = abs((x1 - x2) / sqrt((pow(s1, 2) / n1) + (pow(s2, 2) / n2)))

        #Calculate the degrees of freedom
        num = pow(((pow(s1, 2) / n1) + (pow(s2,2) / n2)), 2)
        den = (pow(pow(s1, 2) / n1, 2) / (n1 - 1)) + (pow(pow(s2, 2) / n2, 2) / (n2 - 1))
        dof = num / den
        
        #Return the two-tailed p value
        return t.sf(t_value, dof) * 2

    def growth_expect_value(self, dict_mutation):
        """Calculate the expectation value for growth selection"""

        #Calculate p(e) or f-1
        prob_expt = (pow(dict_mutation['sel_adj_counts'], (1 / self.growth_gp)) - 1)

        #Prob Sum
        prob_sum = 0
    
        #Loop the rows (gp)
        for i in range (0, int(self.growth_gp) + 1):

            #Get the row list
            row = self.list_pascal[i]

            #Add one
            prob_sum = prob_sum + 1

            #Loop the columns (j = exponent), skip the first column
            for j in range (1, len(row)):

                #Get the pascal number
                pascal_number = self.list_pascal[i][j]
            
                #Add to our running prob sum
                prob_sum = prob_sum + (pascal_number * pow(prob_expt, j))

        return dict_mutation['ref_adj_counts'] * int(prob_sum)

    def fitness_error(self, dict_mutation):
        """Calculate the error based on Poisson noise for the fitness metric"""

        #Select our metric
        if self.metric == "e-wt":
            #e-wt is the variance of the variant added to the variance of the wild-type
            fitness_error = sqrt(pow(dict_mutation['log2error'], 2) + pow(self.dict_summary['log2_wildtype_error'], 2))

        elif self.metric == "growth":
            #Calculate the section for the variant
            section_one = (pow(dict_mutation['log2error'], 2) * 
                            pow((1 / log(2)) * 
                            (1 / (dict_mutation['log2'] + self.growth_gp)), 2))

            #Calculate the section for the wild-type
            section_two = (pow(self.dict_summary['log2_wildtype_error'], 2) * 
                            pow((1 / log(2)) * 
                            (1 / (self.wt_log2 + self.growth_gp)), 2))

            #Combine the sections
            fitness_error = sqrt(section_one + section_two)

        elif self.metric == "facs":
            #Check if the mutant log2 is within the theoretical
            if dict_mutation['log2'] >= self.e_max_theo_facs:
                mut_log2 = self.e_max_theo_facs - 0.001
            else:
                mut_log2 = dict_mutation['log2']

            #FACS is based on flourescence between ei and wt
            section_one = pi * pow(self.facs_pc, 2) * pow(self.facs_sd, 2)

            #Section for the mutation
            section_two = (pow(dict_mutation['log2error'], 2) *
                           pow(pow(2, (mut_log2 + 0.5)) *
                           pow(e, pow(erfinv(1 - self.facs_pc * pow(2, (mut_log2 + 1))), 2)), 2))

            #Section for the wild-type sequence
            section_three = (pow(self.dict_summary['log2_wildtype_error'], 2) *
                             pow(-1 * pow(2, (self.wt_log2 + 0.5)) *
                             pow(e, pow(erfinv(1 - self.facs_pc * pow(2, (self.wt_log2 + 1))), 2)), 2))

            #Combine the sections and report the standard devation
            fitness_error = sqrt(section_one * (section_two + section_three))

        else:
            #Error out
            fitness_error = "NaN"

        return fitness_error

    def wt_synon_sd(self):
        """Process the wild-type synonymous mutation enrichments"""

        #Populate a dict with the counts using the DNA seq as the key
        #0-RefCounts, 1-RefFraction, 2-SelCounts, 3-SelFraction

        #Calculate the log2 enrichment
        list_log2_all = []
        list_log2_12 = []
        list_log2_30 = []
        list_log2_50 = []

        #Loop the dict and calculate the log2 ratio
        for codon in self.dict_wtsynon:
            if self.dict_wtsynon[codon][1] != "NaN" and self.dict_wtsynon[codon][3] != "NaN":
                list_log2_all.append(log(self.dict_wtsynon[codon][3]/self.dict_wtsynon[codon][1], 2))

                #Seperate by read count
                if self.dict_wtsynon[codon][0] >= 50:
                    list_log2_50.append(log(self.dict_wtsynon[codon][3]/self.dict_wtsynon[codon][1], 2))

                if self.dict_wtsynon[codon][0] >= 30:
                    list_log2_30.append(log(self.dict_wtsynon[codon][3]/self.dict_wtsynon[codon][1], 2))

                if self.dict_wtsynon[codon][0] >= 12:
                    list_log2_12.append(log(self.dict_wtsynon[codon][3]/self.dict_wtsynon[codon][1], 2))

        #Calculate our wild-type standard deviation for each metric
        try:
            if self.metric == "growth":
                list_metric_all = [self.fitness_metric_growth(esynon) for esynon in list_log2_all]
                list_metric_12 = [self.fitness_metric_growth(esynon) for esynon in list_log2_12]
                list_metric_30 = [self.fitness_metric_growth(esynon) for esynon in list_log2_30]
                list_metric_50 = [self.fitness_metric_growth(esynon) for esynon in list_log2_50]

            elif self.metric == "facs":

                #Calculate the fitness of the fariants
                list_metric_all = [self.fitness_metric_facs(esynon) for esynon in list_log2_all]
                list_metric_12 = [self.fitness_metric_facs(esynon) for esynon in list_log2_12]
                list_metric_30 = [self.fitness_metric_facs(esynon) for esynon in list_log2_30]
                list_metric_50 = [self.fitness_metric_facs(esynon) for esynon in list_log2_50]

                #Check for inf or -inf and quit
                if float("inf") in list_metric_all or float("-inf") in list_metric_all:
                    print("[Fitness Error] Check your config file for FACS options. Do you have the correct decimal places?")
                    quit()

            else:
                #Calculate for e-wt and use as a placeholder for e/ewt (however wont be used as it is non-gaussian)
                list_metric_all = list_log2_all
                list_metric_12 = list_log2_12
                list_metric_30 = list_log2_30
                list_metric_50 = list_log2_50
        except StatisticsError:
            print("[Fitness] Can not calculate the wild-type synon enrichment standard deviation.")

        #Calculate the standard deviation of the data
        try:
            list_ret_all = [stdev(list_metric_all), len(list_metric_all), mean(list_metric_all), list_metric_all]
        except StatisticsError:
            list_ret_all = [0, 0, 0, []]

        try:
            list_ret_12 = [stdev(list_metric_12), len(list_metric_12), mean(list_metric_12), list_metric_12]
        except StatisticsError:
            list_ret_12 = [0, 0, 0, []]

        try:
            list_ret_30 = [stdev(list_metric_30), len(list_metric_30), mean(list_metric_30), list_metric_30]
        except StatisticsError:
            list_ret_30 = [0, 0, 0, []]

        try:
            list_ret_50 = [stdev(list_metric_50), len(list_metric_50), mean(list_metric_50), list_metric_50]
        except StatisticsError:
            list_ret_50 = [0, 0, 0, []]

        return {'ALL':list_ret_all, '12':list_ret_12, '30':list_ret_30, '50':list_ret_50}

    def output_synon_csv(self, dict_wtsynon_sd):
        """Take the dict of wt_synon values and make a csv report"""

        #Open the output file and write it
        with open(self.out_prefix + "_fitness_synon.csv", "w") as csv_out:

            #Write the header
            csv_out.write("All_Counts,12_Counts,30_Counts,50_Counts\n")

            #Loop and write
            for i in range(0, len(dict_wtsynon_sd['ALL'][3])):
                
                if i >= len(dict_wtsynon_sd['12'][3]):
                    str_12 = ""
                else:
                    str_12 = dict_wtsynon_sd['12'][3][i]

                if i >= len(dict_wtsynon_sd['30'][3]):
                    str_30 = ""
                else:
                    str_30 = dict_wtsynon_sd['30'][3][i]

                if i >= len(dict_wtsynon_sd['50'][3]):
                    str_50 = ""
                else:
                    str_50 = dict_wtsynon_sd['50'][3][i]

                #Write to the file
                csv_out.write(','.join(
                    map(str, [dict_wtsynon_sd['ALL'][3][i], str_12, str_30, str_50, "\n"])))

        return "[Fitness] Fitness values of synon WT codons CSV written"

    def single_ssm_backcalculate(self):
        """Calculate a smaller dict just for SSM single mutations from the larger mutation dict."""

        #Define our mutation types
        mutation_types = '*FWYPMILVAGCSTNQDEHKR'

        #Pre-fill an empty dict with NNK data for all mutated sites
        #{Location : { Mutation { "":X,}}}
        dict_ssm = {}
        for group in self.list_mutation_design:
            for location in group:
                dict_ssm[location] = {} #Create a empty dict per location

                for mut_type in mutation_types:
                    dict_ssm[location][mut_type] = {} #Create a empty dict per location-mutation

                    #If WT add in the WT info else add in the nonsynon data
                    if self.wtaa[location-1] == mut_type:
                        dict_ssm[location][mut_type] = {
                            "location":location,
                            "wt_residue":self.wtaa[location-1],
                            "mutation":mut_type,
                            "reference_counts":self.dict_summary['ref_wt_counts'],
                            "selected_counts":self.dict_summary['sel_wt_counts'],
                            "adjusted_ref_counts":self.dict_summary['ref_wt_counts'],
                            "adjusted_sel_counts":self.dict_summary['sel_wt_counts'],
                            "reference_fraction":self.dict_summary['ref_wtfraction'],
                            "selected_fraction":self.dict_summary['sel_wtfraction'],
                            "log2_enrichment":self.wt_log2,
                            "enrichment_error":self.dict_summary['log2_wildtype_error'],
                            "fitness":0.00,
                            "fitness_error":"NaN",
                            "sd_from_wt":0.00,
                            "welch_p_val":"NaN",
                            }
                    else:
                        #Nonsynonymous mutations
                        dict_ssm[location][mut_type] = {
                            "location":location,
                            "wt_residue":self.wtaa[location-1],
                            "mutation":mut_type,
                            "reference_counts":"0",
                            "selected_counts":"0",
                            "adjusted_ref_counts":"0",
                            "adjusted_sel_counts":"0",
                            "reference_fraction":"NaN",
                            "selected_fraction":"NaN",
                            "log2_enrichment":"NaN",
                            "enrichment_error":"NaN",
                            "fitness":"NaN",
                            "fitness_error":"NaN",
                            "sd_from_wt":"NaN",
                            "welch_p_val":"NaN",
                            }

        #Fill in our data from the self.dict_accepted dict
        for mut_key in self.dict_accepted:
            
            #Back calculate the mutation and location from self.dict_accepted[2] (currently as P,E,E,P...)
            #The location starts at zero       
            adj_location = list(map(int, self.dict_accepted[mut_key]['location'].split(',')))
            adj_mutation = self.dict_accepted[mut_key]['mutation'].split(',')

            #Iterate from the first element to the last element (using the numbering from the dict)
            #Cast as a list with one number and then slice it [0], ==> Returns an integer
            diff_loc = [i for i in range(0, len(adj_location))
                        if adj_mutation[i] != self.wtaa[adj_location[0] + i]][0]           

            #Mutations that have counts in both populations and above a thereshold
            location = adj_location[0] + diff_loc + 1
            mut_type = adj_mutation[diff_loc]

            dict_ssm[location][mut_type] = {
                "location":location,
                "wt_residue":self.wtaa[adj_location[0] + diff_loc],
                "mutation":mut_type,
                "reference_counts":self.dict_accepted[mut_key]['ref_counts'],
                "selected_counts":self.dict_accepted[mut_key]['sel_counts'],
                "adjusted_ref_counts":self.dict_accepted[mut_key]['ref_adj_counts'],
                "adjusted_sel_counts":self.dict_accepted[mut_key]['sel_adj_counts'],
                "reference_fraction":self.dict_accepted[mut_key]['ref_fraction'],
                "selected_fraction":self.dict_accepted[mut_key]['sel_fraction'],
                "log2_enrichment":self.dict_accepted[mut_key]['log2'],
                "enrichment_error":self.dict_accepted[mut_key]['log2error'],
                "fitness":self.dict_accepted[mut_key]['fitness'],
                "fitness_error":self.dict_accepted[mut_key]['fitness_error'],
                "sd_from_wt":self.dict_accepted[mut_key]['sd_from_wt'],
                "welch_p_val":self.dict_accepted[mut_key]['welch_p_val']
                }

        return dict_ssm

    def output_fitness_tsv_single(self, fit_tsv):
        """Output a tsv for the entire dataset (single point mutations)."""
      
        tsv_output = ""
        #Loop each mutation
        for location in fit_tsv:
            for mutation in fit_tsv[location]:

                #Mutations that have counts in both populations and above a thereshold
                tsv_output = tsv_output + '\t'.join(map(str, [
                location,
                fit_tsv[location][mutation]["wt_residue"],
                mutation,
                fit_tsv[location][mutation]["reference_counts"],
                fit_tsv[location][mutation]["selected_counts"],
                fit_tsv[location][mutation]["adjusted_ref_counts"],
                fit_tsv[location][mutation]["adjusted_sel_counts"],
                fit_tsv[location][mutation]["reference_fraction"],
                fit_tsv[location][mutation]["selected_fraction"],
                fit_tsv[location][mutation]["log2_enrichment"],
                fit_tsv[location][mutation]["enrichment_error"],
                fit_tsv[location][mutation]["fitness"],
                fit_tsv[location][mutation]["fitness_error"],
                fit_tsv[location][mutation]["sd_from_wt"],
                str(fit_tsv[location][mutation]["welch_p_val"])+"\n"]))

        #Output our fits and Ewt, then write the tsv files
        tsv_header = '\t'.join(["Location",
                            "WT_Residue",
                            "Mutation",
                            "Reference_Counts",
                            "Selected_Counts",
                            "Adjusted_Ref_Counts",
                            "Adjusted_Sel_Counts",
                            "Reference_Fraction",
                            "Selected_Fraction",
                            "Log2_Enrichment",
                            "Enrichment_Error",
                            "Fitness",
                            "Fitness_Error",
                            "#SD_From_WT",
                            "Welch_PVal\n"])

        #Write the output tsv file
        with open(self.out_prefix + "_Fitness.tsv", 'w') as outfile:
            outfile.write(tsv_header + tsv_output)

        return "[Fitness] Fitness TSV Written"

    def output_fitness_heat(self, dict_ssm):
        """Output a heat for the entire dataset (single mutations only)."""

        #Build the individual sections of the output csv
        mutation_types = '*FWYPMILVAGCSTNQDEHKR'
        output_numbering = ""
        output_wtresi = ""
        for group in self.list_mutation_design:
            for location in group:
                output_numbering = output_numbering +","+ str(location)
                output_wtresi = output_wtresi +","+ self.wtaa[location-1]

        heat_output = ""
        #Iterate the list of mutations
        for aa in mutation_types:
            heat_output = heat_output + aa + ","
            
            #Iterate the groups
            for group in self.list_mutation_design:             
                
                #Iterate each location in the group
                heat_output = heat_output + ','.join(
                    list(str(dict_ssm[location][aa]['fitness'])
                         if location in dict_ssm and aa in dict_ssm[location] else "NaN"
                         for location in group))

            heat_output = heat_output + "\n"

        #Output a summary file of the individual mutations
        with open(self.out_prefix + '_heatmap.csv', 'w') as outfile:
            outfile.write('\n'.join([output_numbering, output_wtresi, heat_output]))

        return "[Fitness] Fitness Heatmap Written"

    def output_fitness_tsv_multi(self):
        """Output a tsv for the entire dataset (multiple mutations)."""

        tsv_output = ""
        #Loop each mutation
        for mut_key in self.dict_accepted:
            #Convert the python location to the actual protein location, and get a wild-type residue
            
            location_adj = ""
            wildtype_resi = ""

            #Get the wildtype residues and then the corrected location
            for location in self.dict_accepted[mut_key]['location'].split(','):
                wildtype_resi = wildtype_resi + self.wtaa[int(location)]+","
                location_adj = location_adj + str(int(location) + 1)+","

            #Mutations that have counts in both populations and above a thereshold
            tsv_output = tsv_output + '\t'.join(map(str, [           
            self.dict_accepted[mut_key]['design'],
            location_adj.rstrip(","),
            wildtype_resi.rstrip(","),
            self.dict_accepted[mut_key]['mutation'],
            self.dict_accepted[mut_key]['ref_counts'],
            self.dict_accepted[mut_key]['sel_counts'],
            self.dict_accepted[mut_key]['ref_adj_counts'],
            self.dict_accepted[mut_key]['sel_adj_counts'],
            self.dict_accepted[mut_key]['ref_fraction'],
            self.dict_accepted[mut_key]['sel_fraction'],
            self.dict_accepted[mut_key]['log2'],
            self.dict_accepted[mut_key]['log2error'],
            self.dict_accepted[mut_key]['fitness'],
            self.dict_accepted[mut_key]['fitness_error'],
            str(self.dict_accepted[mut_key]['sd_from_wt'])+"\n"]))

        #Output our fits and Ewt, then write the tsv files
        tsv_header = '\t'.join(["Designed_Codon",
                            "Location",
                            "WT_Residues",
                            "Mutations",
                            "Reference_Counts",
                            "Selected_Counts",
                            "Adjusted_Ref_Counts",
                            "Adjusted_Sel_Counts",
                            "Reference_Fraction",
                            "Selected_Fraction",
                            "Log2_Enrichment",
                            "Enrichment_Error",
                            "Fitness",
                            "Fitness Error",
                            "#SD_From_WT\n"])

        #Write the output tsv file
        with open(self.out_prefix + "_Fitness.tsv", 'w') as outfile:
            outfile.write(tsv_header + tsv_output)

        return "[Fitness] Fitness TSV Written"

    def fitness(self):
        """Calculate the fitness metric"""

        #Calculate the standard devation of the WT synon mutations
        dict_wtsynon_sd = self.wt_synon_sd()

        #Parse our synon table to make our output string
        output_string = ''.join(map(str, [
                         "All counts SD: ", round(dict_wtsynon_sd['ALL'][0], 4), 
                         " Mean: ", round(dict_wtsynon_sd['ALL'][2], 4),
                         " N=", dict_wtsynon_sd['ALL'][1], "\n",
                         ">=12 counts SD: ", round(dict_wtsynon_sd['12'][0], 4), 
                         " Mean: ", round(dict_wtsynon_sd['12'][2], 4),
                         " N=", dict_wtsynon_sd['12'][1], "\n",
                         ">=30 counts SD: ", round(dict_wtsynon_sd['30'][0], 4), 
                         " Mean: ", round(dict_wtsynon_sd['30'][2], 4),
                         " N=", dict_wtsynon_sd['30'][1], "\n",
                         ">=50 counts SD: ", round(dict_wtsynon_sd['50'][0], 4), 
                         " Mean: ", round(dict_wtsynon_sd['50'][2], 4),
                         " N=", dict_wtsynon_sd['50'][1], "\n",
                         "Wild-Type Log2 Enrichment ", round(self.wt_log2, 4), "\n",
                         "Metric: " + self.metric + "\n",                     
                         ]))

        #Write our raw WT synon codons to a file
        output_string = output_string + self.output_synon_csv(dict_wtsynon_sd) + "\n"

        #Calculate the fitness metric
        for mut_key in self.dict_accepted:

            #Only run mutations that are there
            if self.dict_accepted[mut_key]['log2'] != "NaN":

                #Select our metric
                if self.metric == "e-wt":
                    #e-wt is the log2 of the mutation minus the wild-type enrichment
                    fitness = self.dict_accepted[mut_key]['log2'] - self.wt_log2

                elif self.metric == "growth":
                    #Growth is based on the number of doublings
                    fitness = self.fitness_metric_growth(self.dict_accepted[mut_key]['log2'])

                elif self.metric == "facs":
                    #FACS is based on flourescence between ei and wt
                    fitness = self.fitness_metric_facs(self.dict_accepted[mut_key]['log2'])

                else:
                    #Error out with a non-fitness type
                    print("[Fitness Error] Unknown fitness mode.")
                    quit()

                #Append the fitness and sd_from_wt
                self.dict_accepted[mut_key]['fitness'] = fitness

                #Calculate the fitness error
                self.dict_accepted[mut_key]['fitness_error'] = self.fitness_error(self.dict_accepted[mut_key])

                #How to handle if we don't have any synonymous mutations
                if dict_wtsynon_sd['ALL'][0] != 0 and fitness != "NaN":
                    try:
                        self.dict_accepted[mut_key]['sd_from_wt'] = fitness / dict_wtsynon_sd['ALL'][0]
                    except TypeError:
                        print("[Fitness Error] Unexpected fitness or wt synon values:")
                        print("[Fitness Error] Fitness Value: " + str(fitness))
                        print("[Fitness Error] WT Synon Value: " + str(dict_wtsynon_sd['ALL'][0]))
                        print("[Fitness Error] sd_from_wt set to NaN")
                        self.dict_accepted[mut_key]['sd_from_wt'] = "NaN"

                    #Calculate the Welch p-value between the synonymous wt mutations and non-synon mutations
                    if self.evalue_type == "facs" and sf_import and erfinv_import:
                        #The expectation value is the total cells through the machine * ref fraction
                        evalue_facs = self.evalue_facs_cellcount * self.dict_accepted[mut_key]['ref_fraction']

                        #Calculate the p-value from the T-Test with Welch's correction
                        self.dict_accepted[mut_key]['welch_p_val'] = self.welch_pval(
                            dict_wtsynon_sd['ALL'][1], evalue_facs,
                            dict_wtsynon_sd['ALL'][0], self.dict_accepted[mut_key]['fitness_error'],
                            0, fitness)

                    elif self.evalue_type == "growth" and sf_import:
                        #The expectation value is based on cell doubling and using pascals triangle
                        evalue_growth = self.growth_expect_value(self.dict_accepted[mut_key])

                        #Calculate the p-value from the T-Test with Welch's correction
                        self.dict_accepted[mut_key]['welch_p_val'] = self.welch_pval(
                            dict_wtsynon_sd['ALL'][1], evalue_growth,
                            dict_wtsynon_sd['ALL'][0], self.dict_accepted[mut_key]['fitness_error'],
                            0, fitness)

                    else:
                        self.dict_accepted[mut_key]['welch_p_val'] = "NaN"

                else:
                    self.dict_accepted[mut_key]['sd_from_wt'] = "NaN"
                    self.dict_accepted[mut_key]['welch_p_val'] = "NaN"

            else:
                #Append NaN to fitness if we don't have it
                self.dict_accepted[mut_key]['fitness'] = "NaN"
                self.dict_accepted[mut_key]['sd_from_wt'] = "NaN"
                self.dict_accepted[mut_key]['welch_p_val'] = "NaN"
                self.dict_accepted[mut_key]['fitness_error'] = "NaN"

        #Output a heatmap or residue frequency depending on mutation type
        if self.library_type == "single":
            #Backcalculate to a more constrained library
            dict_ssm = self.single_ssm_backcalculate()

            #Calculate the residue frequency
            output_string = output_string + self.output_fitness_heat(dict_ssm) + "\n"

            #Output the tsv file
            output_string = output_string + self.output_fitness_tsv_single(dict_ssm) + "\n"

            #Open and write the output SSM file
            
            output_string = output_string + save_pact_file(dict_ssm, self.out_prefix + '_' + "SSM_fitness") + "\n"

        elif self.library_type == "multiple":
            #Output the tsv file
            output_string = output_string + self.output_fitness_tsv_multi() + "\n"

        #Save a .pact file for the wt synon dict
        output_string = output_string + save_pact_file(dict_wtsynon_sd, self.out_prefix + '_' + "fitness_wtsynon") + "\n"

        #Save a .pact file for the fitness dict
        output_string = output_string + save_pact_file(self.dict_accepted, self.out_prefix + '_' + "fitness_nonsynon") + "\n"

        #Print our output
        print(output_string)

        return output_string

if __name__ == '__main__':

    #Remind the user that the module needs to be ran within the context of PACT
    print("[Fitness Error] This module needs to be ran within the context of PACT.")
