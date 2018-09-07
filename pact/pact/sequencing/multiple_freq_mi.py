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

"""Calculate the multiple residue codon frequency and mutual information"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Multiple Residue Freq and MI Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from configparser import NoOptionError
from itertools import chain, product, combinations, groupby
import numpy as np
from pact.pact_common import save_pact_file, open_pact_file

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class multiple_freq_mi:
    """Calculate the multiple residue codon frequency and mutual information"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Check if the config file section is defined
        if not self.config_file.has_section("multiple_freq_mi"):
            print("[Multi Freq and MI Error] the config file is missing the section [multiple_freq_mi]")
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

        self.mutation_types = '*FWYPMILVAGCSTNQDEHKR'
        
        #Import the library mode (Single or Multiple)
        if self.dict_protocolconfig['library_type'].lower() != 'multiple':
            print("[Multi Freq and MI Error] Unknown library type (expect multiple).")
            quit()

        #Import settings from the config file on file imports
        try:
            #Test if we have a special input file
            if len(self.config_file.get('multiple_freq_mi', 'pact_fitness_nonsynon')) > 0:
                file_pact_accepted = self.config_file.get('multiple_freq_mi', 'pact_fitness_nonsynon')
            else:
                file_pact_accepted = self.out_prefix + "_fitness_nonsynon.pact"

            #Check our inported files and load
            self.dict_accepted = open_pact_file(file_pact_accepted.rstrip('.pact'))

        except NoOptionError:
            print("[Multi Freq and MI Error] The fitness config file is incorrect.")
            print("[Multi Freq and MI Error] Missing an option flag that starts with file_.")
            quit()

        #Should we only accept mutations with log2 values?
        try:
            if self.config_file.get('multiple_freq_mi', 'frequency_log2_filter').lower() == "true":
                self.log2_requirement = True
            else:
                self.log2_requirement = False
        except:
            self.log2_requirement = False

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
                    print("[Multi Freq and MI Error] The codon list is missing either a start or end point around the 'n' location.")
                    quit()
            else:
                #Add the other ranges without changing them
                self.list_mutation_design.append(group)

        #Flatten the list of mutation designs, sort it, and remove dupes
        self.flatten_locations = set(sorted([int(loc) for group in self.list_mutation_design for loc in group]))

        return

    def pactdict_to_nparray(self):
        """Convert a PACT dict to a dict with np arrays"""

        #Numpy arrays do not like to be appended to, 
        #therefore regular lists in Python need to be used.

        #Setup our dicts
        dict_ref = {}
        dict_sel = {}

        #Loop the pact dict
        #Note the key is a string (it looks like a list)
        for key in self.dict_accepted:

            #Check for the design group (matches the config file, count from 1)
            if self.dict_accepted[key]['design'] not in dict_ref:
                dict_ref[self.dict_accepted[key]['design']] = []
                dict_sel[self.dict_accepted[key]['design']] = []

            #Reject mutations that pass our filter (i.e. lack a log2) if requested
            #Will have to think if this is needed?? Saves a if statement
            #if self.log2_requirement and self.dict_accepted[key]['log2'] == "NaN":
            #    continue

            #Append the mutation to our list, multiply by the number of counts
            dict_ref[self.dict_accepted[key]['design']].append([self.dict_accepted[key]['mutation']] 
                                                               * self.dict_accepted[key]['ref_counts'])
            dict_sel[self.dict_accepted[key]['design']].append([self.dict_accepted[key]['mutation']] 
                                                               * self.dict_accepted[key]['sel_counts'])

        #Convert to numpy array, filtering out empty, flatten list
        # The list so far will have a empty if read count = 0
        # The list will also have mutations like [['muts', 'muts', 'muts', 'muts'], ['muts']...
        # from_iterable will return a list that removes emptys and flattens the entries
        # Result will be ['muts', 'muts', 'muts', 'muts']
        # Then will have to split on the comma to make [['m,u,t,s'], ['m,u,t,s']...
        # That will be our numpy array
        for key in dict_ref:

            #Append our np array
            dict_ref[key] = np.array([muts.split(',') 
                                              for muts in list(chain.from_iterable(dict_ref[key]))])
            dict_sel[key] = np.array([muts.split(',') 
                                              for muts in list(chain.from_iterable(dict_sel[key]))])

        return dict_ref, dict_sel

    def nparray_process(self, dict_nparray_ref, dict_nparray_sel):
        """Calculate the per residue mutation frequency, log2 enrichment"""

        #Setup a dict to work into
        dict_sitewise = {}

        #Loop the design locations    
        for key in dict_nparray_ref:

            #Get the length of the design location (the key is a string)
            list_sites = key.lstrip('[').rstrip(']').split(',')

            #Transpose the array (convert cols to rows)

            ##FIX: Check if we want to keep tolist (converts numpy to list???)
            ref_transpose = np.array(dict_nparray_ref[key]).T.tolist()
            sel_transpose = np.array(dict_nparray_sel[key]).T.tolist()

            #Get the total count
            total_unique_ref, total_counts_ref = np.unique(ref_transpose, return_counts=True)
            total_unique_sel, total_counts_sel = np.unique(sel_transpose, return_counts=True)
            ref_total_sum = np.sum(total_counts_ref)
            sel_total_sum = np.sum(total_counts_sel)

            #Loop each column (i = site)
            for i in range(0, len(list_sites)):

                #Location
                location = int(list_sites[i])

                #If the site has alread been added
                if location in dict_sitewise:
                    print("[Multiple Freq and MI Error] Site " + list_sites[i] +
                         " has already been added, is this site present in more than one design group?")
                else:
                    dict_sitewise[location] = {}

                #Get the counts
                unique_ref, counts_ref = np.unique(ref_transpose[i], return_counts=True)
                unique_sel, counts_sel = np.unique(sel_transpose[i], return_counts=True)

                #Get the column sum
                ref_column_sum = np.sum(counts_ref)
                sel_column_sum = np.sum(counts_sel)

                #Get the site-wise frequency
                ref_site_freq = np.divide(counts_ref, ref_column_sum)
                sel_site_freq = np.divide(counts_sel, sel_column_sum)

                #Get the change in frequency
                freq_change = np.subtract(sel_site_freq, ref_site_freq)

                #Get the group-wise frequency
                ref_group_freq = np.divide(counts_ref, ref_total_sum)
                sel_group_freq = np.divide(counts_sel, sel_total_sum)
                
                #Get the log2 ratio
                log2 = np.log2(np.divide(sel_group_freq, ref_group_freq))

                #Add our information to our dict (AminoAcid, 
                #dict_sitewise[list_sites[i]] = np.asarray((unique_ref, counts_ref, counts_sel,
                #                                          ref_site_freq, sel_site_freq, freq_change,
                #                                          ref_group_freq, sel_group_freq, log2)).T

                #Add our information
                for j in range(0, len(unique_ref)):

                    #Add to our dict
                    dict_sitewise[location][unique_ref[j]] = {
                        'location':location,
                        "wt_residue":self.wtaa[location - 1],
                        "mutation":unique_ref[j],
                        "reference_counts":counts_ref[j],
                        "selected_counts":counts_sel[j],                        
                        "ref_site_freq":ref_site_freq[j],
                        "sel_site_freq":sel_site_freq[j],
                        "site_freq_change":freq_change[j],
                        "reference_fraction":ref_group_freq[j],
                        "selected_fraction":sel_group_freq[j],
                        "log2_enrichment":log2[j],                          
                        }

                #Normalize to the wild-type residue at the site
                for j in range(0, len(unique_ref)):
                    
                    dict_sitewise[location][unique_ref[j]]['fitness'] = (
                        dict_sitewise[location][unique_ref[j]]["log2_enrichment"] - 
                        dict_sitewise[location][self.wtaa[location - 1]]["log2_enrichment"]
                        )

        return dict_sitewise

    def mutual_information(self, dict_nparray_sel):
        """Calculate the mutual information between residues"""

        #The input is a dict with the key of the sites and an array of sequences
        #Programming note: nparr_filtered == dict_nparray_sel[design_key]

        #Loop the designs, design_key = list of sites
        for design_key in dict_nparray_sel:

            print("[Mutual Information] Mutual information for design group:")
            print(design_key)

            #Locations (convert str to list)
            list_locations = literal_eval(design_key)

            #Store a list of numpy arrays (index = site)(0 = letter, 1 = frequency)
            list_single_freq = []

            #Get the single residue freq
            for j in range(0, len(list_locations)):
            
                #Count the unique aminos
                unique, counts = np.unique(dict_nparray_sel[design_key].T[j], return_counts=True)

                #Get the site-wise frequency
                site_freq = np.divide(counts, np.sum(counts))

                #Append to a #3D list 1D = everything, 2D = location, 3D = [mutation, freq]
                list_single_freq.append(np.asarray((unique, site_freq)).T.tolist())

            #Create a list of tuples (#, #) that have combinations of sites choose 2
            list_sites = list(combinations(range(0, len(list_locations)), 2))

            #Return a 2D list
            list_combinations = []

            #Site locations, tuple = (0, 1)
            for mapping in list_sites:

                #Get the single mutation frequencies
                nparr_single_freq1 = np.array(list_single_freq[mapping[0]])
                nparr_single_freq2 = np.array(list_single_freq[mapping[1]])

                #Get the list of column variants from the dict
                #Slice on all rows only return the column we want
                #Reshape into a 2D list
                nparr_col1 = dict_nparray_sel[design_key][:,mapping[0]].reshape((dict_nparray_sel[design_key].shape[0], 1))
                nparr_col2 = dict_nparray_sel[design_key][:,mapping[1]].reshape((dict_nparray_sel[design_key].shape[0], 1))

                #Return an array of the double variants, transpose and slice to get 2d from 3d array
                nparr_dbl_variants = np.asarray((nparr_col1, nparr_col2)).T[0]

                #Return the unique double pairings and the counts
                nparr_dbl_combos, nparr_dbl_counts = np.unique(nparr_dbl_variants, axis=0, return_counts=True)

                #Create the single freq lists
                #Replace with a numpy solution in the future
                #https://stackoverflow.com/questions/29320873/replace-values-of-a-numpy-array-by-values-from-another-numpy-array
                list_freq1 = []
                list_freq2 = []
                for combo in nparr_dbl_combos:
                    list_freq1.append(nparr_single_freq1[np.where(nparr_single_freq2[:,0] == combo[0])][0][1])
                    list_freq2.append(nparr_single_freq2[np.where(nparr_single_freq2[:,0] == combo[1])][0][1])

                #Add to our lists
                #Replace with a numpy solution in the future
                for i in range(0, nparr_dbl_combos.shape[0]):
                    list_combinations.append([
                        mapping[0],
                        mapping[1],
                        nparr_dbl_combos[i][0],
                        nparr_dbl_combos[i][1],
                        list_freq1[i],
                        list_freq2[i],
                        nparr_dbl_counts.reshape((nparr_dbl_combos.shape[0], 1))[i][0],
                        dict_nparray_sel[design_key].shape[0]
                        ])

            #Convert to numpy array
            nparr_combinations = np.array(list_combinations)  

            #Extract the frequency and counts, reshape to 2D, and type
            freq_one = np.asarray(nparr_combinations.T[4].reshape((nparr_combinations.shape[0], 1))).astype(float)
            freq_two = np.asarray(nparr_combinations.T[5].reshape((nparr_combinations.shape[0], 1))).astype(float)
            counts_doubles = np.asarray(nparr_combinations.T[6].reshape((nparr_combinations.shape[0], 1))).astype(int)
            counts_total = np.asarray(nparr_combinations.T[7].reshape((nparr_combinations.shape[0], 1))).astype(int)
            
            #Calculate the frequency of the double variant
            frequency_doubles = np.divide(counts_doubles, counts_total)

            #Calculate the MI
            mi = np.log2(np.divide(frequency_doubles, np.multiply(freq_one, freq_two)))

            #Calculate the APC (Average Product Correction) Dunn Bioinformatics
            #APC = Avg MI of site n to all sites except to itself, Avg MI of site m to all except itself
            #Divided by the mean MI of all sites (except sites to itself)
            #Subtract this from the MI

            #Iterate each row of location combinations
            list_apc = []
            for row in nparr_combinations[:,0:2].tolist():
                #Use isin twice, the first to see if the row is in the 2d array
                #Then a second time to see if True is in the 2d array (the first
                #will return a array of [True False] or other combo
                nparr_mi_avg_idx = np.any(np.isin(nparr_combinations[:,0:2], row), axis=1)

                #Calculate the APC value
                list_apc.append(np.divide(np.multiply(mi[nparr_mi_avg_idx].mean(), mi[nparr_mi_avg_idx].mean()), mi.mean()))
            
            #Calculate the MI-APC per site    
            nparr_mi_apc = np.subtract(mi, np.array(list_apc).reshape((nparr_combinations.shape[0], 1)))
           
            #Output a file
            str_header = "Site1,Site2,Resi1,Resi2,Site1Freq,Site2Freq,Counts,TotalCounts,MI-APC"
            str_filename = "_" + str(list_locations[0]) + "_" + str(list_locations[-1])
            np.savetxt(self.out_prefix + str_filename + "_MI.csv", 
                       np.hstack((nparr_combinations, nparr_mi_apc)), 
                       fmt='%s,%s,%s,%s,%s,%s,%s,%s,%s', 
                       header=str_header)
        return

    def output_heat(self, dict_sitewise):
        """Output a heatmap of calculated information"""     

        #Build the individual sections of the output csv
        output_numbering = ""
        output_wtresi = ""
        for location in self.flatten_locations:
            output_numbering = output_numbering +","+ str(location)
            output_wtresi = output_wtresi +","+ self.wtaa[location - 1]

        #Build the main sections
        output_log2 = ""
        output_delta_freq = ""

        #Loop the amino acids
        for aa in self.mutation_types:

            #Write the header
            output_log2 = output_log2 + aa + ","
            output_delta_freq = output_delta_freq + aa + ","

            #Loop the locations
            for location in self.flatten_locations:

                output_log2 = output_log2 + str(dict_sitewise[location][aa]['log2_enrichment']) + ","
                output_delta_freq = output_delta_freq + str(dict_sitewise[location][aa]['site_freq_change']) + ","

            #Write the ending of the line
            output_log2 = output_log2 + "\n"
            output_delta_freq = output_delta_freq + "\n"

        #Output a summary file of the individual mutations
        with open(self.out_prefix + '_SingleResidueFreqs.csv', 'w') as outfile:
            outfile.write('\n'.join(["Log2 Enrichment",
                                     output_numbering,
                                     output_wtresi,
                                     output_log2,
                                     "Change in Frequency",
                                     output_numbering,
                                     output_wtresi,
                                     output_delta_freq]))

        return "[Multi Freq and MI] Wrote: " + self.out_prefix + '_SingleResidueFreqs.csv'

    def output_tsv(self, dict_sitewise):
        """Output a tsv of calculated information"""
      
        tsv_output = ""
        #Loop each mutation
        for location in self.flatten_locations:
            for mutation in self.mutation_types:

                #Mutations that have counts in both populations and above a thereshold
                tsv_output = tsv_output + ','.join(map(str, [
                location,
                dict_sitewise[location][mutation]["wt_residue"],
                mutation,
                dict_sitewise[location][mutation]["reference_counts"],
                dict_sitewise[location][mutation]["selected_counts"],
                dict_sitewise[location][mutation]["ref_site_freq"],
                dict_sitewise[location][mutation]["sel_site_freq"],
                dict_sitewise[location][mutation]["site_freq_change"],
                dict_sitewise[location][mutation]["reference_fraction"],
                dict_sitewise[location][mutation]["selected_fraction"],
                dict_sitewise[location][mutation]["log2_enrichment"],
                str(dict_sitewise[location][mutation]["fitness"])+"\n"]))

        #Output our fits and Ewt, then write the tsv files
        tsv_header = ','.join(["Location",
                            "WT_Residue",
                            "Mutation",
                            "Reference_Counts",
                            "Selected_Counts",
                            "Ref_Site_Freq",
                            "Sel_Site_Freq",
                            "Site_Freq_Change",
                            "Ref_Group_Freq",
                            "Sel_Group_Freq",
                            "Log2_Enrichment",
                            "Fitness\n"])

        #Write the output tsv file
        with open(self.out_prefix + "_Multi_Freq.csv", 'w') as outfile:
            outfile.write(tsv_header + tsv_output)

        return "[Multi Freq and MI] Fitness TSV Written"

    def multiple_freq_mi(self):
        """Main routine to handle multi-site libraries"""

        #Convert a PACT dict to a np array (returns a dict with design as key and array of ['m', 'u', 't', 's'], ['m'..
        print("[Multiple Freq and MI] Converting the PACT file to a numpy array.")
        dict_nparray_ref, dict_nparray_sel = self.pactdict_to_nparray()

        #Calculate the residue frequency
        print("[Multiple Freq and MI] Processing the numpy array for frequency and log2.")
        dict_sitewise = self.nparray_process(dict_nparray_ref, dict_nparray_sel)

        #Output our count, freq, and log2 information
        print("[Multiple Freq and MI] Writing heatmap, column, and ssm.pact.")
        print(self.output_heat(dict_sitewise))
        self.output_tsv(dict_sitewise)
        print(save_pact_file(dict_sitewise, self.out_prefix + '_multiple_SSM'))

        #Calculate the mutual information (considering residues within groups)
        print("[Multiple Freq and MI] Mutual Information.")
        self.mutual_information(dict_nparray_sel)

        #Print our output
        output_string = "Wrote per-residue CSV\n"
        #print(output_string)

        return output_string

if __name__ == '__main__':

    #Remind the user that the module needs to be ran within the context of PACT
    print("[Multi Freq and MI Error] This module needs to be ran within the context of PACT.")
