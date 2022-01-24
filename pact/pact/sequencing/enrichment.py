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

"""enrichment - calculate the enrichment of mutations"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Enrichment Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from math import log, sqrt, exp, pow
from pact.pact_common import file_checker

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class enrichment:
    """Calculate the enrichment of mutations"""

    def __init__(self, settings_dict):
        """Initialize the class varibles"""
        
        #Assign general DNA and AA settings
        self.wtdna = settings_dict['WTDNA']
        self.FirstAAMutated = int(settings_dict['FirstAAMutated'])-1
        self.LastAAMutated = int(settings_dict['LastAAMutated'])-1
        self.WTDNARegion = self.wtdna[self.FirstAAMutated*3:(self.LastAAMutated + 1) * 3]

        #Setup the output file prefix
        self.out_prefix = settings_dict['Out_Prefix']

        #Check our inported files
        if file_checker(settings_dict["Ref_Count_WildType"]):
            self.file_ref_count_wt = settings_dict["Ref_Count_WildType"]
        
        if file_checker(settings_dict["Sel_Count_WildType"]):
            self.file_sel_count_wt = settings_dict["Sel_Count_WildType"]
        
        if file_checker(settings_dict["Ref_Count"]):
            self.file_ref_count = settings_dict["Ref_Count"]        
            
        if file_checker(settings_dict["Sel_Count"]):
            self.file_sel_count = settings_dict["Sel_Count"]           
            
        if file_checker(settings_dict["Ref_Count_Rejected"]):
            self.file_ref_count_rejected = settings_dict["Ref_Count_Rejected"]        
            
        if file_checker(settings_dict["Sel_Count_Rejected"]):
            self.file_sel_count_rejected = settings_dict["Sel_Count_Rejected"]         

        #Set our thresholds
        try:
            self.ref_count_threshold = int(settings_dict["Ref_Count_Threshold"])
            self.sel_count_threshold = int(settings_dict["Sel_Count_Threshold"])
        except ValueError:
            print("[Error] The read count threshold is incorrectly set, defaulting to 12.")
            self.ref_count_threshold = 12
            self.sel_count_threshold = 12

        #Do we enforce a strict count threshold?
        try:
            if settings_dict["Strict_Count_Threshold"].lower() == "true":
                self.strict_threshold = True
            else:
                self.strict_threshold = False
        except:
            self.strict_threshold = False

        #
        # do we consider mutations rejected in our design yet pass
        # our read count filters to be in the total library count value?
        # 2022.1 - default is now false as theoretically a large amount
        # of non-designed mutations could skew the dataset distribution
        # in the non-normalized enrichments
        #
        try:
            if settings_dict["consider_rejected"].lower() == "true":
                self.consider_rejected = True
            else:
                self.consider_rejected = False
        except:
            self.consider_rejected = False
            
        #Import the mutation design list
        self.list_mutation_design = []
        for group in literal_eval(settings_dict['mutcodons']):
            #If we have n between two numbers
            if 'n' in group:
                #We need exactly three to expand the range
                if len(group) == 3:
                    #Get the start and end position, then iterate through the two points
                    self.list_mutation_design.append(list(i for i in range(group[0], group[-1] + 1)))
                else:
                    print("[Enrichment Error] The codon list is missing either a start or end point around the 'n' location.")
                    quit()
            else:
                #Add the other ranges without changing them
                self.list_mutation_design.append(group)
        return

    def wt_synon_filter(self, list_lines):
        """Filter the wild-type reads for mutations within our codons (called by wt_synon)"""

        list_synon = []

        #Go through each line
        for line in list_lines: 
            #Split the line
            splitline = line.split('\t')

            #Zip the DNA sequence and translated sequence
            joined_dna = list(zip(splitline[1], self.WTDNARegion))

            #Count the number of mismatches from the input DNA
            #End here if we are at the primary sequence
            if sum(1 for a, b in joined_dna if a != b) == 0:
                continue

            #Get the codon locations of the mismatches
            list_mutated_codons = []
            for i in range(0, len(joined_dna)):
                if joined_dna[i][0] != joined_dna[i][1]: #If the DNA base doesn't match...
                    list_mutated_codons.append(divmod(i, 3)[0] + self.FirstAAMutated) #Find the codon number using divmod

            #Now check to see if the mutation is within a designed codon
            list_designmatch = []
            for codon in list_mutated_codons: #Loop through the number of mismatches
                list_designmatch.append(None) #Add a none for each codon

                for j in range(0, len(self.list_mutation_design)): #Loop through the number of designed sets
                    for k in self.list_mutation_design[j]: #Get each designed AA number
                        if codon == k:
                            list_designmatch.pop() #Remove a none from the end
                            list_designmatch.append(self.list_mutation_design[j]) #Add the proper location

            #Now lets filter the matches
            if list_designmatch.count(list_designmatch[0]) == len(list_designmatch) and list_designmatch[0] is not None:
                list_synon.append([splitline[1], int(splitline[2])])                    

        return list_synon

    def wt_synon(self, list_reference, list_selected, counts_reftotal, counts_seltotal):
        """Calculate the WT synon enrichment"""

        #Map the ref and sel WT synon lists
        dict_wtsynon_enrichment = {}

        #Populate a dict with the counts using the DNA seq as the key
        #0-RefCounts, 1-RefFraction, 2-SelCounts, 3-SelFraction
        for codon in self.wt_synon_filter(list_reference):
            dict_wtsynon_enrichment[codon[0]] = [codon[1], None, None, None]

        for codon in self.wt_synon_filter(list_selected):
            try:
                dict_wtsynon_enrichment[codon[0]][2] = codon[1]
            except KeyError:
                dict_wtsynon_enrichment[codon[0]] = [None, None, codon[1], None]

        #Calculate the fraction
        for codon in dict_wtsynon_enrichment:
            #Calculate the reference
            if dict_wtsynon_enrichment[codon][0] is not None:
                dict_wtsynon_enrichment[codon][1] = dict_wtsynon_enrichment[codon][0]/counts_reftotal
            else:
                dict_wtsynon_enrichment[codon][1] = "NaN"

            #Calculate the selected
            if dict_wtsynon_enrichment[codon][2] is not None:
                dict_wtsynon_enrichment[codon][3] = dict_wtsynon_enrichment[codon][2]/counts_seltotal
            else:
                dict_wtsynon_enrichment[codon][3] = "NaN"
       
        return dict_wtsynon_enrichment

    def count_importer(self, file_input, list_columns):
        """Read in the counted mutation file into a dict and count the total number of mutations."""

        #Open the file with the counted mutations        
        with open(file_input, 'r') as infile:
            list_lines = infile.read().splitlines()[1:]

        #Read the file into a dict
        dict_mutations = {}
        for line in list_lines:
            splitline = line.split("\t")
            design = splitline[list_columns[0]]

            #Check if we have a header
            if design == "ID":
                continue

            location = splitline[list_columns[1]]
            mutation = splitline[list_columns[2]]
            counts = int(splitline[list_columns[3]])

            #Check for existing values
            if location+mutation in dict_mutations:
                dict_mutations[location+mutation][3] = dict_mutations[location+mutation][3] + counts
            else:
                dict_mutations[location+mutation] = [design, location, mutation, counts]

        return dict_mutations, list_lines

    def dict_merge_count(self, dict_reference, dict_selected):
        """Merge the reference and selected dicts."""

        #Prepare the dict with 0's for selected counts (add column #4)
        for mut_key in dict_reference:
            dict_reference[mut_key].append(0)

        #Add in selected counts
        for mut_key in dict_selected:
            #Check if the key exists in the reference dict
            if mut_key in dict_reference:
                dict_reference[mut_key][4] = dict_selected[mut_key][3]
            else:
                dict_reference[mut_key] = [dict_selected[mut_key][0], dict_selected[mut_key][1],
                                        dict_selected[mut_key][2], 0, dict_selected[mut_key][3]]

        #Add in virtual counts (for when a variant fell out or blew up and it's partner doesn't exist)
        count_ref = 0
        count_sel = 0
        for mut_key in dict_reference:
            #Handle if we are enforcing strict read counts
            if self.strict_threshold:           
                #Add in 0 if below our threshold
                if dict_reference[mut_key][3] < self.ref_count_threshold:
                    dict_reference[mut_key].append(0)
            
                #Else add in the original value
                else:
                    dict_reference[mut_key].append(dict_reference[mut_key][3])

          
                #Add in 0 if both are below our threshold (column #6) 
                if dict_reference[mut_key][4] < self.sel_count_threshold:
                    dict_reference[mut_key].append(0)
            
                #Else add in the original value
                else:
                    dict_reference[mut_key].append(dict_reference[mut_key][4])

            else:
                #Add in 1 for the reference if the selected is above our threshold (column #5)
                if dict_reference[mut_key][4] >= self.sel_count_threshold and dict_reference[mut_key][3] == 0:
                    dict_reference[mut_key].append(1)
            
                #Add in 0 if both are below our threshold
                elif dict_reference[mut_key][4] < self.sel_count_threshold and dict_reference[mut_key][3] < self.ref_count_threshold:
                    dict_reference[mut_key].append(0)
            
                #Else add in the original value
                else:
                    dict_reference[mut_key].append(dict_reference[mut_key][3])


                #Add in 1 for the selected if the reference is above our threshold (column #6)
                if dict_reference[mut_key][3] >= self.ref_count_threshold and dict_reference[mut_key][4] == 0:
                    dict_reference[mut_key].append(1)
            
                #Add in 0 if both are below our threshold
                elif dict_reference[mut_key][4] < self.sel_count_threshold and dict_reference[mut_key][3] < self.ref_count_threshold:
                    dict_reference[mut_key].append(0)
            
                #Else add in the original value
                else:
                    dict_reference[mut_key].append(dict_reference[mut_key][4])

            #Get the total number of counts above a threshold
            count_ref = count_ref + dict_reference[mut_key][5]
            count_sel = count_sel + dict_reference[mut_key][6]

        return dict_reference, count_ref, count_sel

    def fraction_log2(self, dict_selection, counts_reftotal, counts_seltotal):
        """Calculate the count fraction and log2 ratio"""
        #0design, 1location, 2mutation, 3ref_counts, 4sel_counts, 5ref_adj_counts, 6sel_adj_counts
        #7ref_fraction, 8sel_fraction, 9log2, 10log2error

        #Calculate the fraction
        for mut_key in dict_selection:
            #Add a new value 7:ref_fraction
            if dict_selection[mut_key][5] != 0:
                dict_selection[mut_key].append(dict_selection[mut_key][5]/counts_reftotal)
            else:
                dict_selection[mut_key].append("NaN")

            #Add a new value 8:sel_fraction
            if dict_selection[mut_key][6] != 0:
                dict_selection[mut_key].append(dict_selection[mut_key][6]/counts_seltotal)
            else:
                dict_selection[mut_key].append("NaN")
        
            #Add a new value 9:log2
            if dict_selection[mut_key][7] != "NaN" and dict_selection[mut_key][8] != "NaN":
                dict_selection[mut_key].append(log(dict_selection[mut_key][8] / dict_selection[mut_key][7], 2))
            else:
                dict_selection[mut_key].append("NaN")

            #Add a new value 10:log2 error 1SD
            if dict_selection[mut_key][5] != 0 and dict_selection[mut_key][6] != 0:
                dict_selection[mut_key].append(sqrt(pow(log(exp(1), 2), 2)*((1/dict_selection[mut_key][5]) + (1/dict_selection[mut_key][6]))))
            else:
                dict_selection[mut_key].append("NaN")

        return dict_selection

    def save_result_dict(self, dict_save, dict_name, convert_to_dict=False):
        """Save the results of this program into a python friendly form."""

        #Import the pickle dump method
        from pickle import dump, HIGHEST_PROTOCOL

        #Do we convert the dict of lists to dict of dicts?
        if convert_to_dict:
            #Since the dict mappings are the same for the nonsynon, reject, and fitness we can remap them
            #0design, 1location, 2mutation, 3ref_counts, 4sel_counts, 5ref_adj_counts, 6sel_adj_counts
            #7ref_fraction, 8sel_fraction, 9log2, 10log2error
            #Append the fitness as list item #11 FM and #12 Num of WT SD away

            new_dict_save = {}

            for entry in dict_save:
                #Add a new dict for the entry
                new_dict_save[entry] = {}

                #Check to see if the fitness values exist
                dict_to_add = {
                    'design':dict_save[entry][0],
                    'location':dict_save[entry][1],
                    'mutation':dict_save[entry][2],
                    'ref_counts':dict_save[entry][3],
                    'sel_counts':dict_save[entry][4],
                    'ref_adj_counts':dict_save[entry][5],
                    'sel_adj_counts':dict_save[entry][6],
                    'ref_fraction':dict_save[entry][7],
                    'sel_fraction':dict_save[entry][8],
                    'log2':dict_save[entry][9],
                    'log2error':dict_save[entry][10]
                    }
                    
                #Add the new dict
                new_dict_save[entry] = dict_to_add

            #Open and write the output file
            with open(self.out_prefix + '_' + dict_name + '.pact', 'wb') as outfile:
                dump(new_dict_save, outfile, protocol=HIGHEST_PROTOCOL)
        else: 
            #Open and write the output file
            with open(self.out_prefix + '_' + dict_name + '.pact', 'wb') as outfile:
                dump(dict_save, outfile, protocol=HIGHEST_PROTOCOL)

        return '[Enrichment] Output dumped as ' + self.out_prefix + '_' + dict_name + '.pact'

    def enrichment(self):
        """Calculate the mutation enrichment"""

        #The column number from the counted file for design, location, mutation, and counts
        list_rejected_cols = [0, 1, 2, 3] #List columns for rejected: [0:design = 0, 1:location = 1, 2:mutation = 2, 3:counts = 3]
        list_nonsynon_cols = [0, 1, 3, 8] #List columns for nonsynon: [0:design = 0, 1:location = 1, 2:mutation = 3, 3:counts = 8]
        list_wildtype_cols = [0, 0, 0, 2] #List columns for wildtype: [0:design = 0, 1:location = 0, 2:mutation = 0, 3:counts = 2]

        # Read in the counted mutation files
        # 2022.1 update for rejected mutations
        if self.consider_rejected:
            dict_reject_ref, list_reject_ref = self.count_importer(self.file_ref_count_rejected, list_rejected_cols)
            dict_reject_sel, list_reject_sel = self.count_importer(self.file_sel_count_rejected, list_rejected_cols)
        
        dict_wildtype_ref, list_wt_ref = self.count_importer(self.file_ref_count_wt, list_wildtype_cols)
        dict_wildtype_sel, list_wt_sel = self.count_importer(self.file_sel_count_wt, list_wildtype_cols)
        
        dict_nonsynon_ref, list_nonsynon_ref = self.count_importer(self.file_ref_count, list_nonsynon_cols)
        dict_nonsynon_sel, list_nonsynon_sel = self.count_importer(self.file_sel_count, list_nonsynon_cols)

        # merge the rejected and nonsynon selections into one dict each
        # 2022.1 update for rejected mutations
        if self.consider_rejected:
            dict_reject, counts_rejected_ref, counts_rejected_sel = self.dict_merge_count(dict_reject_ref, dict_reject_sel)
        
        dict_nonsynon, counts_nonsynon_ref, counts_nonsynon_sel = self.dict_merge_count(dict_nonsynon_ref, dict_nonsynon_sel)

        #Count the wildtype
        counts_wildtype_ref = sum(dict_wildtype_ref[count][3] for count in dict_wildtype_ref)
        counts_wildtype_sel = sum(dict_wildtype_sel[count][3] for count in dict_wildtype_sel)

        #
        # 2022.1 update - updated handling of strict counts
        # As read depth increases, the rejected counts could skew the non-normalized dataset
        # Therefore, it is suggested to not consider and remove the rejected counts
        # from the total read count, (to retain the old default set a value of True for consider_rejected
        # in the [enrichment] config section
        # 
        if self.consider_rejected:
            # sum the counts
            counts_reftotal = counts_rejected_ref + counts_wildtype_ref + counts_nonsynon_ref
            counts_seltotal = counts_rejected_sel + counts_wildtype_sel + counts_nonsynon_sel
        else:
            # default action
            
            # sum the counts without rejected mutations adding to the total read amount
            counts_reftotal = counts_wildtype_ref + counts_nonsynon_ref
            counts_seltotal = counts_wildtype_sel + counts_nonsynon_sel
            
            # set rejected counts to zero
            counts_rejected_ref = 0
            counts_rejected_sel = 0

        #Calculate the wild-type enrichment
        try:
            fraction_refwt = counts_wildtype_ref / counts_reftotal
            fraction_selwt = counts_wildtype_sel / counts_seltotal
            log2_wildtype = log(fraction_selwt / fraction_refwt, 2)
            log2wt_error = sqrt(pow(log(exp(1), 2), 2)*((1/counts_wildtype_ref) + (1/counts_wildtype_sel)))
            str_wt_status = "[Enrichment] Calculated the wild-type log2 enrichment.\n"
        except ValueError:
            str_wt_status = ("[Error] Math domain error in calculation of the WT enrichment value.\n"
                  "[Error] Did you use the same fastq file for the fwd or rev directions?\n"
                  "[Error] Is wild-type not present?\n"
                  "[Error] Log2 WT and frequency is set to 'Error' in the dataset.\n"
                  )
            fraction_refwt = "Error"
            fraction_selwt = "Error"
            log2_wildtype = "Error"
            log2wt_error = "Error"

        #Get the synonymous mutation SD
        dict_wtsynon_enrichment = self.wt_synon(list_wt_ref, list_wt_sel, counts_reftotal, counts_seltotal)

        #[0design, 1location, 2mutation, 3ref_counts, 4sel_counts, 5ref_adj_counts, 6sel_adj_counts]

        # calculate the fractions and log2 enrichment
        # 2022.1 update for rejected mutations
        if self.consider_rejected:
            dict_reject = self.fraction_log2(dict_reject, counts_reftotal, counts_seltotal)
        
        dict_nonsynon = self.fraction_log2(dict_nonsynon, counts_reftotal, counts_seltotal)

        #0design, 1location, 2mutation, 3ref_counts, 4sel_counts, 5ref_adj_counts, 6sel_adj_counts
        #7ref_fraction, 8sel_fraction, 9log2, 10log2error

        #Assemble a dict that is a summary of the most important information
        dict_enrichment = {
            'log2_wildtype':log2_wildtype,
            'log2_wildtype_error':log2wt_error,
            'ref_wtfraction':fraction_refwt,
            'sel_wtfraction':fraction_selwt,
            'ref_wt_counts':counts_wildtype_ref,
            'sel_wt_counts':counts_wildtype_sel,
            'ref_accepted_counts':counts_nonsynon_ref,
            'sel_accepted_counts':counts_nonsynon_sel,
            'ref_reject_counts':counts_rejected_ref,
            'sel_reject_counts':counts_rejected_sel,
            'ref_total_counts':counts_reftotal,
            'sel_total_counts':counts_seltotal,
            'ref_count_threshold':self.ref_count_threshold,
            'sel_count_threshold':self.sel_count_threshold,
            }

        #Assemble an output string
        output_string = ''.join(map(str, [
            "[Enrichment] Processed the wild-type, accepted, and rejected reads.\n",
            "[Enrichment] Calculated the log2 enrichments of all above threshold mutations.\n",
            str_wt_status,
            "[Enrichment] Wild-Type Log2 Enrichment ", log2_wildtype, "\n",
            "[Enrichment] Counts Ref Wild-Type Synonymous ", counts_wildtype_ref, "\n",
            "[Enrichment] Counts Sel Wild-Type Synonymous ", counts_wildtype_sel, "\n",
            "[Enrichment] Counts Ref Accepted Non-Synonymous ", counts_nonsynon_ref, "\n",
            "[Enrichment] Counts Sel Accepted Non-Synonymous ", counts_nonsynon_sel, "\n",
            "[Enrichment] Counts Ref Rejected Non-Synonymous ", counts_rejected_ref, "\n",
            "[Enrichment] Counts Sel Rejected Non-Synonymous ", counts_rejected_sel, "\n",
            "[Enrichment] Total Reference Counts ", counts_reftotal, "\n",
            "[Enrichment] Total Selected Counts ", counts_seltotal, "\n"]))

        #Save our calculated enrichment summary (counts, wt log2, etc)
        output_string = output_string + self.save_result_dict(dict_enrichment, "enrichment_summary") + "\n"

        #Save our synonymous dict
        output_string = output_string + self.save_result_dict(dict_wtsynon_enrichment, "enrichment_wtsynon") + "\n"

        # Save our rejected dict
        # 2022.1 update for rejected mutations
        if self.consider_rejected:
            output_string = output_string + self.save_result_dict(dict_reject, "enrichment_reject_nonsynon", True) + "\n"
        
        #Save our accepted dict
        output_string = output_string + self.save_result_dict(dict_nonsynon, "enrichment_accept_nonsynon", True) + "\n"

        #Output our string
        print(output_string)

        return output_string

if __name__ == '__main__':

    #Remind the user that the module needs to be ran within the context of PACT
    print("[Enrichment Error] This module needs to be ran within the context of PACT.")
