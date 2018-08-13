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

"""library_stats - sequencing statistics of the library"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Library Stats Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from pickle import load
from operator import mul
from functools import reduce
from statistics import stdev, mean, median, StatisticsError
from pact.pact_common import file_checker

#Import the codon_analysis method from pact.tools.codon_condenser
try:
    from pact.tools.codon_condenser import codon_condenser
except ImportError:
    print("[Library Stats Error] pact.tools.codon_condenser was not found, exiting.")
    quit()

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class library_stats:
    """Calculate the library stats."""

    def __init__(self, settings_dict):
        """Import our settings and saved config file."""

        #Assign general DNA and AA settings
        self.wtdna = settings_dict['WTDNA']
        self.wtaa = settings_dict['WTAA']
        self.FirstAAMutated = int(settings_dict['FirstAAMutated'])-1
        self.LastAAMutated = int(settings_dict['LastAAMutated'])-1
        self.out_prefix = settings_dict['Out_Prefix']

        #Import the library mode (Single or Multiple)
        if settings_dict['library_type'].lower() == 'single':
            self.library_type = "single"
        elif settings_dict['library_type'].lower() == 'multiple':
            self.library_type = "multiple"
        else:
            print("[Library Stats Error] Unknown mode (expect single or multiple).")
            quit()

        #Check our inported files and load
        if file_checker(settings_dict["file_summary"]):
            with open(settings_dict['file_summary'], 'rb') as file_summary:
                self.dict_summary = load(file_summary)

        if file_checker(settings_dict["pact_fitness_nonsynon"]):
            with open(settings_dict['pact_fitness_nonsynon'], 'rb') as file_accepted:
                self.fitness = load(file_accepted)        

        if file_checker(settings_dict["pact_fitness_wtsynon"]):
            with open(settings_dict['pact_fitness_wtsynon'], 'rb') as file_wtsynon:
                self.dict_wtsynon = load(file_wtsynon)

        #Parse our pact files for global variables
        self.ref_threshold = int(self.dict_summary['ref_count_threshold'])
        self.sel_threshold = int(self.dict_summary['sel_count_threshold'])

        #Let the mutation design list
        self.list_mutation_design = []
        for group in literal_eval(settings_dict['mutcodons']):
            #If we have n between two numbers
            if 'n' in group:
                #We need exactly three to expand the range
                if len(group) == 3:
                    #Get the start and end position, then iterate through the two points
                    self.list_mutation_design.append(list(i for i in range(group[0], group[-1] + 1)))
                else:
                    print("[Error] The codon list is missing either a start or end point around the 'n' location.")
                    quit()
            else:
                #Add the other ranges without changing them
                self.list_mutation_design.append(group)

        #'codon_type':"{'NNK':[[161,'n',256]]}"
        #Let the mutation design list
        list_codon_design = []
        dict_codon_design = literal_eval(settings_dict['codon_type'])
        
        #Enumerate each codon key
        for codon in dict_codon_design:

            #Enumerate the list
            for group in dict_codon_design[codon]:

                #If we have n between two numbers
                if 'n' in group:

                    #We need exactly three to expand the range
                    if len(group) == 3:
                        #Get the start and end position, then iterate through the two points
                        list_codon_design.append(
                            {i: codon for i in range(group[0], group[-1] + 1)})
                    else:
                        print("[Error] The codon list is missing either a start or end point around the 'n' location.")
                        quit()
                else:
                    #Add the other ranges without changing them
                    list_codon_design.append({j: codon for j in group})

        #Now we need to flatten and combine the dicts
        self.dict_codon_design = {k: v for list_entry in list_codon_design for k, v in list_entry.items()}

        return
 
    def import_codon_design(self, codon):
        """For a given designed codon, return which amino acids are encoded."""

        #Create our object
        obj_cc = codon_condenser()

        #Return the list of codons applicable for our design
        #{'A': ['GCC'], 'C': ['TGC'], 'D': ['GAC'], ...
        dict_codons, output_str = obj_cc.codon_analysis(codon)

        #Make a list of the amino acids that have codons present
        list_aminos = list(key for key, value in dict_codons.items() if len(dict_codons[key]) > 0)

        #Return the designed codons and the list of aminos expected
        return dict_codons, list_aminos
   
    def theo_coverage(self):
        """Calculate the theoretical amino acid coverage of the entire experiment."""

        #Loop the input design per region
        #[[161, 162, 163], [167, 169, 171], [175, 177, 179],...

        #Get the number of amino acids at each codon (keeping the grouping structure)
        #[[21, 21, 21], [21, 21, 21], ...
        #If the user has an incorrect codon_type varible set I.E. doesn't match the tile number then a keyerror will show here

        try:
            list_amino_counts = list([len(self.import_codon_design(self.dict_codon_design[codon_num])[1])
                                    for codon_num in region]
                                    for region in self.list_mutation_design)
        except KeyError:
            print("[Library Stats Error] Double check the codon_type variable as a residue number is incorrect.")
            print("[Library Stats Error] Please verify the config file and re-run. Quitting on failure.")
            quit()

        #Get the specific amino acids at each codon in each group
        list_codon_aminos = list([self.import_codon_design(self.dict_codon_design[codon_num])[1]
                                for codon_num in region]
                                for region in self.list_mutation_design)

        #Get the specific codons at each group
        list_codon_codons = list([self.import_codon_design(self.dict_codon_design[codon_num])[0]
                                for codon_num in region]
                                for region in self.list_mutation_design)

        #Get the total number of codons at each position
        list_codon_counts = list([sum([len(codon[amino]) for amino in codon])
                                for codon in group]
                                for group in list_codon_codons)

        #For SSM do addition, for multiple get the product for different grouping stats
        if self.library_type == "single":
            #Treat each codon as independent without combinations with other codons
            #Returns the sum of each group
            list_group_aminos = list(sum(region) for region in list_amino_counts)

            #Returns the sum of each group (gets the total number of codons)
            list_group_codons = list(sum(region) for region in list_codon_counts)

        else:
            #Treat each codon in combinations with other codons in each design group
            #Return the product of each group (reduce function copys the mul function which is x*y over the list)
            list_group_aminos = list(reduce(mul, region) for region in list_amino_counts)

            #Returns the sum of each group (gets the total number of codons)
            list_group_codons = list(reduce(mul, region) for region in list_codon_counts)

        #Exports the combined total
        count_total_aminos = sum(list_group_aminos) #Total amino acid counts
        count_total_codons = sum(list_group_codons) #Total codon counts

        #Zip the lists together [Design, Group Count]
        list_zip_design_group = list(zip(self.list_mutation_design, list_group_aminos))

        return (count_total_aminos, count_total_codons, list_group_aminos,
                count_total_aminos, list_zip_design_group, list_codon_aminos)

    def nonsynon_coverage(self):
        """Calculate the amino acid coverage of the library for a given design."""

        #Count the total number of accepted nonsynon mutations (with a log2 enrichment)
        count_nonsyon = sum(1 for entry in self.fitness if self.fitness[entry]['log2'] != "NaN")
        
        return count_nonsyon

    def nonsense(self):
        """Return the average stop codon enrichment"""

        #Get a list of enrichments of variants with a stop codon
        list_nonsense = list(self.fitness[entry]['log2'] for entry in self.fitness
                            if self.fitness[entry]['log2'] != "NaN" and self.fitness[entry]['mutation'].find('*') != -1)

        #Return the mean and SD of stop codons enrichment
        try:
            return mean(list_nonsense), stdev(list_nonsense)
        except StatisticsError:
            print("[Library Stats Error] Not enough reads, setting the Median and Mean of stop codons to 0.")
            return ("NaN", "NaN")

    def read_count_hist(self, list_counts, count_median, count_mean, threshold, str_filename):
        """Plot the dist of the read counts"""

        #Import the graphing toolkit
        try:
            from matplotlib import use, rcParams
        except ImportError:
            print("[Library Stats Error] matplotlib is not installed, will not plot read distributions.")
            return

        #Set the matplotlib base settings (headless, font)
        use('Agg')
        rcParams['font.sans-serif']='Arial'

        from matplotlib import pyplot as plt
        from numpy import logspace, ones_like

        #Make the x-axis bins and y-axis weights (to make it out of 1)
        bins = logspace(0, 5, 100) #Make a list with 100 bins between 1e0 and 1e5
        weights = ones_like(list_counts)/float(len(list_counts)) #Weight the y-axis values to make it out of 1

        #Create the histogram
        fig, ax = plt.subplots()
        counts, bins, patches = ax.hist(list_counts, bins, weights=weights, color='black')

        #Label the graph
        plt.xlabel("Number of reads per variant")
        plt.ylabel("Fraction of total reads")
    
        #Add the scale, mean, and median lines
        plt.xscale('log')
        plt.axvline(count_mean, ymin=0, ymax=1, color='blue', linestyle='--')
        plt.axvline(count_median, ymin=0, ymax=1, color='red', linestyle='--')
        plt.axvline(threshold, ymin=0, ymax=1, color='orange', linestyle='--')
    
        #Add the legend
        textstr = "Mean: "+str(count_mean)+"\nMedian: "+str(count_median)+"\nN: "+str(len(list_counts))
        plt.text(0.8, 0.95, textstr, transform=ax.transAxes, fontsize=12, 
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))

        #Save the figure as a png
        plt.savefig(self.out_prefix + str_filename + "_count_hist.png", dpi=600)

        return

    def library_stats(self):
        """Produce a readout of the total library stats."""

        #Get the read counts into a list
        list_refadj_counts = list(self.fitness[mut]['ref_adj_counts'] for mut in self.fitness 
                                  if self.fitness[mut]['ref_adj_counts'] != 0)
        list_seladj_counts = list(self.fitness[mut]['sel_adj_counts'] for mut in self.fitness 
                                  if self.fitness[mut]['sel_adj_counts'] != 0)

        #Calculate the median and mean read counts for the adjusted populations
        try:
            median_refadj_readcount = round(median(list_refadj_counts))
            mean_refadj_readcount = round(mean(list_refadj_counts))
            median_seladj_readcount = round(median(list_seladj_counts))
            mean_seladj_readcount = round(mean(list_seladj_counts))
        except StatisticsError:
            print("[Libary Stats Error] Median or Mean with empty data")
            quit()

        #Create the read count histograms
        self.read_count_hist(list_refadj_counts, median_refadj_readcount, mean_refadj_readcount, self.ref_threshold, "_RefAdj")
        self.read_count_hist(list_seladj_counts, median_seladj_readcount, mean_seladj_readcount, self.sel_threshold, "_SelAdj")

        #Calculate the ref and sel if both are above the threshold
        list_ref_counts = list(self.fitness[mut]['ref_counts'] for mut in self.fitness
                               if self.fitness[mut]['ref_counts'] >= self.ref_threshold 
                               and self.fitness[mut]['sel_counts'] >= self.sel_threshold)
        list_sel_counts = list(self.fitness[mut]['sel_counts'] for mut in self.fitness
                               if self.fitness[mut]['ref_counts'] >= self.ref_threshold 
                               and self.fitness[mut]['sel_counts'] >= self.sel_threshold)

        #Calculate the median and mean read counts for the adjusted populations
        string_output = ""
        try:
            median_ref_readcount = round(median(list_ref_counts))
            mean_ref_readcount = round(mean(list_ref_counts))
            median_sel_readcount = round(median(list_sel_counts))
            mean_sel_readcount = round(mean(list_sel_counts))

            #Create the read count histograms
            self.read_count_hist(list_ref_counts, median_ref_readcount, mean_ref_readcount, self.ref_threshold, "_Ref")
            self.read_count_hist(list_sel_counts, median_sel_readcount, mean_sel_readcount, self.sel_threshold, "_Sel")

        except StatisticsError:
            print("[Library Stats Error] Not enough reads, setting the Median and Mean readcounts to 0.\n")
            string_output = "[Library Stats Error] Not enough reads, setting the Median and Mean readcounts to 0.\n"
            median_ref_readcount = "insufficient reads to perform statistics"
            mean_ref_readcount = "insufficient reads to perform statistics"
            median_sel_readcount = "insufficient reads to perform statistics"
            mean_sel_readcount = "insufficient reads to perform statistics"      

        #Return the theoretical amino acid/codon coverage
        tuple_theo_coverage = self.theo_coverage()

        #Get the total nonsynon coverage with a fitness metric
        count_nonsyon = self.nonsynon_coverage()

        #Get the average and SD of nonsense mutations enrichment
        tuple_nonsense = self.nonsense()

        if tuple_nonsense[0] == "NaN":
            str_nonsense = "insufficient reads to perform statistics"
        else:
            str_nonsense = ("{0:.4f}".format(float(tuple_nonsense[0])) + 
            " +/- " + 
            "{0:.4f}".format(float(tuple_nonsense[1])))

        #Calculate the number of locations in the protein
        int_number_locations = len([item for sublist in self.list_mutation_design for item in sublist])

        #Construct the output string
        string_output = string_output + '\n'.join([
            "[Library Stats] Information from the deep sequencing library:",
            "***Tile Information***",
            "Total Number of Codons: " + str(self.LastAAMutated - self.FirstAAMutated + 1),
            "First Codon: " + str(self.FirstAAMutated + 1) + self.wtaa[self.FirstAAMutated],
            "Last Codon: " + str(self.LastAAMutated + 1) + self.wtaa[self.LastAAMutated],
            "Count threshold for reference: " + str(self.ref_threshold),
            "Count threshold for selected: " + str(self.sel_threshold),
            
            ("Wild-Type Enrichment and SD: " 
             "{0:.4f}".format(float(self.dict_summary['log2_wildtype'])) +
             " +/- " + 
             "{0:.4f}".format(float(self.dict_wtsynon['ALL'][0]))),
            "Extended information for standard deviation of synonymous WT mutations:",

            "All counts SD: " + str(round(self.dict_wtsynon['ALL'][0], 4)), 
            "Mean: " + str(round(self.dict_wtsynon['ALL'][2], 4)),
            "N: " + str(self.dict_wtsynon['ALL'][1]),
            ">=12 counts SD: " + str(round(self.dict_wtsynon['12'][0], 4)), 
            "Mean: " + str(round(self.dict_wtsynon['12'][2], 4)),
            "N: " + str(self.dict_wtsynon['12'][1]),
            ">=30 counts SD: " + str(round(self.dict_wtsynon['30'][0], 4)), 
            "Mean: " + str(round(self.dict_wtsynon['30'][2], 4)),
            "N: " + str(self.dict_wtsynon['30'][1]),
            ">=50 counts SD: " + str(round(self.dict_wtsynon['50'][0], 4)), 
            "Mean: " + str(round(self.dict_wtsynon['50'][2], 4)),
            "N: " + str(self.dict_wtsynon['50'][1]),

            "Average nonsense (stop) enrichment and SD: " + str_nonsense,

            "",
            "***Counts for Reference***",
            "Total Number of Reads: " + str(self.dict_summary['ref_total_counts']) + " (100%)",
            ("Number (%) with No Nonsynonymous Mutations: " + 
             str(self.dict_summary['ref_wt_counts']) + 
             " ({0:.1f}%)".format((self.dict_summary['ref_wt_counts']/self.dict_summary['ref_total_counts'])*100)),
            ("Number (%) with Accepted Nonsynonymous Mutations: " + 
             str(self.dict_summary['ref_accepted_counts']) + 
             " ({0:.1f}%)".format((self.dict_summary['ref_accepted_counts']/self.dict_summary['ref_total_counts'])*100)),
            ("Number (%) with Rejected Nonsynonymous Mutations: " + 
             str(self.dict_summary['ref_reject_counts']) + 
             " ({0:.1f}%)".format((self.dict_summary['ref_reject_counts']/self.dict_summary['ref_total_counts'])*100)),
            "",
            "***Counts for Selected***",
            "Total Number of Reads: " + str(self.dict_summary['sel_total_counts']) + " (100%)",            
            ("Number (%) with No Nonsynonymous Mutations: " + 
             str(self.dict_summary['sel_wt_counts']) + 
             " ({0:.1f}%)".format((self.dict_summary['sel_wt_counts']/self.dict_summary['sel_total_counts'])*100)),
            ("Number (%) with Accepted Nonsynonymous Mutations: " + 
             str(self.dict_summary['sel_accepted_counts']) + 
             " ({0:.1f}%)".format((self.dict_summary['sel_accepted_counts']/self.dict_summary['sel_total_counts'])*100)),
            ("Number (%) with Rejected Nonsynonymous Mutations: " + 
             str(self.dict_summary['sel_reject_counts']) + 
             " ({0:.1f}%)".format((self.dict_summary['sel_reject_counts']/self.dict_summary['sel_total_counts'])*100)),
            "",
            "***Amino Acid Coverage***",
            "Theoretical size of library (WT+Nonsynonymous mutations): " + str(tuple_theo_coverage[0]),
            "Theoretical size of library (Nonsynonymous mutations): " + str(
                tuple_theo_coverage[0] - int_number_locations),
            ("Coverage of nonsynonymous mutations with log2 values: " + 
             str(count_nonsyon) + 
             " ({0:.1f}%)".format(
                 (count_nonsyon / (tuple_theo_coverage[0] - int_number_locations)) * 100)),
            "",
            "***Codon Coverage***",
            "Median Nonsynonymous Read Count of Reference Library (adjusted): " + str(median_refadj_readcount),
            "Median Nonsynonymous Read Count of Selected Library (adjusted): " + str(median_seladj_readcount),
            "Number of nonsynonymous variants (adjusted): " + str(len(list_refadj_counts)),
            "Median Nonsynonymous Read Count of Reference Library (strict threshold): " + str(median_ref_readcount),
            "Median Nonsynonymous Read Count of Selected Library (strict threshold): " + str(median_sel_readcount),
            "Number of nonsynonymous variants (strict threshold): " + str(len(list_ref_counts)),
            "Theoretical size of library (WT+Nonsynonymous codons): " + str(tuple_theo_coverage[1]),

            ("Fold oversampling of total reference reads of codon combinations: " + 
             "{0:.1f}x".format(((self.dict_summary['ref_total_counts']-self.dict_summary['ref_reject_counts']) / tuple_theo_coverage[1]))),

            ("Fold oversampling of total selected reads of codon combinations: " + 
             "{0:.1f}x".format(((self.dict_summary['sel_total_counts']-self.dict_summary['sel_reject_counts']) / tuple_theo_coverage[1]))),

            ("Fold oversampling of Nonsynonymous reference reads of codon combinations: " + 
             "{0:.1f}x".format(((self.dict_summary['ref_accepted_counts'])
                                / (tuple_theo_coverage[1] - int_number_locations)))),

            ("Fold oversampling of Nonsynonymous selected reads of codon combinations: " + 
             "{0:.1f}x".format(((self.dict_summary['sel_accepted_counts']) 
                                / (tuple_theo_coverage[1]-len(self.list_mutation_design[0])))))
            ])

        #Save our output as a file
        with open(self.out_prefix + "_stats.txt", "w") as file_statsout:
            file_statsout.write(string_output)
        
        #Output our workup
        print(string_output)

        return string_output

if __name__ == '__main__':
    #Remind the user that the module needs to be ran within the context of PACT
    print("[Library Stats Error] This module needs to be ran within the context of PACT.")
