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

"""Classification Analysis - Classify mutations and perform a broad sequence and structure analysis used for model development"""

from sys import version_info

#Setup our protocol name
str_protocol_name = "Classification Analysis"

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
from pact.pact_common import save_pact_file, open_pact_file
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class classification_analysis:
    """Classify mutations and perform a broad sequence and structure analysis used for model development."""

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

    def protocol(self):
        """Main entrypoint for the protocol"""

        #Create a output log file that we can append to
        with open(self.directory + self.output_prefix + "_" + strftime("%m_%d_%Y") + "-" +
                 strftime("%H_%M_%S") + '_output.txt', 'w') as file_output:
            file_output.write(self.pact_preamble + "\n")

            """
            *****************************************
            Pact Combine
            *****************************************
            """
            if self.dict_workflow['combinepact']:
                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('combinepact'):           
                    print("[Protocols:" + str_protocol_name + " Error] The combinepact config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [combinepact] section.")
                    quit()

                #Import our combinepact class
                try:
                    from pact.analysis.combine_pact import combine_pact
                except ImportError:
                    print("[Protocols:" + str_protocol_name + " Error] combine_pact was not found.")

                #Create the object then call the merger
                obj_combine = combine_pact(self.obj_cfgparser, self.dict_programs, {})

                #The dict will be like {'dataset name': {data...
                dict_merged_datasets = obj_combine.combine_pact()

            """
            *****************************************
            Classify our mutations
            *****************************************
            """
            #Build a dict_classified with [location][mutation] = "DEL/NEU/BEN/NONE"
            
            #Get the config file elements
            try:
                class_column = self.obj_cfgparser.get("variant_classification", "class_column").lower()
                class_threshold = float(self.obj_cfgparser.get("variant_classification", "class_threshold"))
            except NoOptionError:
                print("[Protocols:" + str_protocol_name + " Error] Missing [variant_classification] config file elements.")
                quit()
            except ValueError:
                print("[Protocols:" + str_protocol_name + " Error] Incorrect [variant_classification] config file elements.")
                quit()
            except TypeError:
                print("[Protocols:" + str_protocol_name + " Error] Incorrect [variant_classification] config file elements.")
                quit()

            #Make a dict to add our classifications into
            dict_classified = {}

            #Classify each dataset
            for dataset in dict_merged_datasets:

                #Add if not existing
                if dataset not in dict_classified:
                    dict_classified[dataset] = {}

                #Loop the locations
                for loc in dict_merged_datasets[dataset]:
                
                    #Add a new location if not in the dict
                    if loc not in dict_classified[dataset]:
                        dict_classified[dataset][loc] = {}

                    #Loop the muts
                    for mut in dict_merged_datasets[dataset][loc]:

                        #Skip WT, stop, and NaN
                        if (mut == dict_merged_datasets[dataset][loc][mut]['wt_residue'] or
                            mut == "*" or
                            dict_merged_datasets[dataset][loc][mut][class_column] == "NaN"):

                            dict_classified[dataset][loc][mut] = "UNCLASSIFIED"
                            continue

                        #Get the fitness value from the dataset
                        mut_value = float(dict_merged_datasets[dataset][loc][mut][class_column])

                        #Assign a classification of deleterious, slightly deleterious, or neutral
                        if mut_value <= (-1 * class_threshold):
                            dict_classified[dataset][loc][mut] = "DEL"
                        elif (mut_value > (-1 * class_threshold) and mut_value < class_threshold):
                            dict_classified[dataset][loc][mut] = "NEU"
                        elif mut_value >= class_threshold:
                            dict_classified[dataset][loc][mut] = "BEN"

                        #if mut_value < -1:
                        #    dict_classified[dataset][loc][mut] = "DEL"
                        #elif mut_value >= -1 and mut_value < -0.3:
                        #    dict_classified[dataset][loc][mut] = "NEU"
                        #elif mut_value >= -0.3:
                        #    dict_classified[dataset][loc][mut] = "BEN"

            """
            *****************************************
            Count the basal classifiers
            *****************************************
            """
            if self.dict_workflow['basal_count']:

                #Import our class
                try:
                    from pact.analysis.basal_count import basal_count
                except ImportError:
                    print("[Protocols:" + str_protocol_name + " Error] pact.analysis.basal_count was not found.")
                    quit()
        
                #Create our object
                obj_basal = basal_count(self.obj_cfgparser, self.dict_programs, {})

                #Count our basal rates
                for dataset in dict_classified:
                    print("[Protocols:" + str_protocol_name + "] Basal Fitness Counts for dataset: " + dataset)
                    file_output.write("[Protocols:" + str_protocol_name + "] Basal Fitness Counts for dataset: " + dataset)

                    file_output.write(obj_basal.basal_count(dict_classified[dataset]) + "\n")

            """
            *****************************************
            DNA Filtering/Alignment or PSSM Object
            *****************************************
            """
            if (self.dict_workflow['blastp_align_filter'] or 
                self.dict_workflow['pssm'] or 
                self.dict_workflow['pssm_reader']):

                #Import Check Output
                from subprocess import check_output

                #Import our class
                try:
                    from pact.analysis.sequence.homology_pssm import homology_classifier
                except ImportError:
                    print("[Protocols:" + str_protocol_name + " Error] pact.analysis.sequence.homology_pssm was not found.")
                    quit()
        
                #Create our object
                obj_homology = homology_classifier(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #Count our classifiers
                print("[Protocols:" + str_protocol_name + "] Homology PSSM")
                file_output.write("[Protocols:" + str_protocol_name + "] Homology PSSM")

            """
            *****************************************
            DNA Filtering/Alignment
            *****************************************
            """
            if self.dict_workflow['blastp_align_filter']:
               
                #Convert our XML file
                print("[Protocols:" + str_protocol_name + "] xml_to_fasta")
                file_output.write("[Protocols:" + str_protocol_name + "] xml_to_fasta\n")
                obj_homology.xml_to_fasta()

                #Run CD-HIT on our new fasta file
                print("[Protocols:" + str_protocol_name + "] cdhit")
                file_output.write("[Protocols:" + str_protocol_name + "] cdhit\n")

                #Check to see if the number of processes is logical
                self.processes = self.obj_cfgparser.get("blastp_align_filter", "processes")
                if int(self.processes) <= 0:
                    self.processes = "2"

                check_output([self.dict_programs['cdhit'],
                              "-i",
                              self.directory + self.output_prefix + ".fa",
                              "-o",
                              self.directory + self.output_prefix + ".afa",
                              "-c",
                              str(self.obj_cfgparser.get("blastp_align_filter", "cdhit_clustering_threshold")),
                              "-M",
                              "40000",
                              "-T",
                              str(self.processes)])

                #Check to see if we have WT in our cdhit output
                print("[Protocols:" + str_protocol_name + "] cdhit_wtcheck")
                file_output.write("[Protocols:" + str_protocol_name + "] cdhit_wtcheck\n")
                obj_homology.cdhit_wt_check()

                #Run MUSCLE on our new fasta file
                print("[Protocols:" + str_protocol_name + "] muscle")
                file_output.write("[Protocols:" + str_protocol_name + "] muscle\n")
                check_output([self.dict_programs['muscle'],
                              "-in",
                              self.directory + self.output_prefix + ".afa",
                              "-out",
                              self.directory + self.output_prefix + ".msa"])

                #Process our MSA (needs to be on for PSIBlast)
                print("[Protocols:" + str_protocol_name + "] processmsa")
                file_output.write("[Protocols:" + str_protocol_name + "] processmsa\n")
                list_msa = obj_homology.process_msa()

                #Save our list
                print("[Protocols:" + str_protocol_name + "] Saving our MSA")
                file_output.write("[Protocols:" + str_protocol_name + "] Saving our MSA\n")
                save_pact_file(list_msa, self.directory + self.output_prefix + '_' + "list_msa")

            """
            *****************************************
            PSSM
            *****************************************
            """
            if self.dict_workflow['pssm']:
                #Open our list
                print("[Protocols:" + str_protocol_name + "] Opening our MSA")
                file_output.write("[Protocols:" + str_protocol_name + "] Opening our MSA\n")
                list_msa = open_pact_file(self.directory + self.output_prefix + '_' + "list_msa")

                #Split our msa for PSIBlast (needs to be on for PSIBlast)
                print("[Protocols:" + str_protocol_name + "] msa_split")
                file_output.write("[Protocols:" + str_protocol_name + "] msa_split\n")
                list_pbcmds = obj_homology.msa_split(list_msa)

                #Run PSIBlast
                print("[Protocols:" + str_protocol_name + "] psiblast")
                file_output.write("[Protocols:" + str_protocol_name + "] psiblast\n")
                for command in list_pbcmds:
                    check_output([self.dict_programs['psiblast'], *command])

                #Import our PSSM data
                print("[Protocols:" + str_protocol_name + "] pssm_file_import")
                file_output.write("[Protocols:" + str_protocol_name + "] pssm_file_import\n")
                dict_pssm = obj_homology.pssm_file_import()

                #Save our heatmap
                print("[Protocols:" + str_protocol_name + "] Saving a PSSM .csv heatmap")
                file_output.write(obj_homology.pssm_output_heat(dict_pssm) + "\n")

                #Save our csv
                print("[Protocols:" + str_protocol_name + "] Saving a PSSM .csv column data")
                file_output.write(obj_homology.pssm_output_csv(dict_pssm) + "\n")

                #Save our PACT File
                print("[Protocols:" + str_protocol_name + "] Saving a PSSM .pact file")
                file_output.write(save_pact_file(dict_pssm, self.directory + self.output_prefix + '_' + "PSSM") + "\n")

            """
            *****************************************
            Read stored PSSM files
            *****************************************
            """
            if (self.dict_workflow['pssm_reader'] or self.dict_workflow['consensus']):
                #Open our PACT File
                print("[Protocols:" + str_protocol_name + "] Opening a PSSM .pact file")
                dict_pssm = open_pact_file(self.directory + self.output_prefix + '_' + "PSSM")

            """
            *****************************************
            PDB Import Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['pdb_import']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('pdb_import'):           
                    print("[Protocols:" + str_protocol_name + " Error] The pdb_import config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [pdb_import] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.pdb_import import pdb_import
                except ImportError:
                    print("[Protocols:" + str_protocol_name + " Error] pdb_import was not found.")

                #Create the object then call the merger
                obj_pdb = pdb_import(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #The dict will be like {'pdb name': {data...
                dict_pdb = obj_pdb.pdb_import()

            """
            *****************************************
            Back to Consensus Analyses
            *****************************************
            """
            if self.dict_workflow['consensus']:
                #Import our class
                try:
                    from pact.analysis.sequence.consensus import consensus
                except ImportError:
                    print("[Protocols:" + str_protocol_name + "] pact.analysis.basal_count was not found.")
                    quit()
        
                #Create our object
                obj_consensus = consensus(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #Get the wild-type sequence information
                dict_wtcons = obj_consensus.wt_consensus(dict_pssm)

            """
            *****************************************
            Residue Chemical/Size
            *****************************************
            """
            if self.dict_workflow['residue_chemical_size']:

                #Import our residue_chemical_size class
                try:
                    from pact.analysis.sequence.residue_chemical_size import residue_chemical_size
                except ImportError:
                    print("[Protocols:Enzyme Solubility Error] residue_chemical_size was not found.")

                #Create the object then call the merger
                obj_rcs = residue_chemical_size(self.obj_cfgparser, self.dict_programs, {})

            """
            *****************************************
            Distance to Active Site
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['distance_to_active']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('distance_to_active'):           
                    print("[Protocols:Enzyme Solubility Error] The distance_to_active config file is incorrect.")
                    print("[Protocols:Enzyme Solubility Error] There is something wrong with the [distance_to_active] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.dist_to_active import dist_to_active
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.dist_to_active was not found.")
                    quit()
        
                #Create our object
                obj_dtoa = dist_to_active(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_dtoa_dist = obj_dtoa.dta_dist(dict_pdb)

            """
            *****************************************
            Contact Number
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['contact_number']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('contact_number'):           
                    print("[Protocols:Enzyme Solubility Error] The contact_number config file is incorrect.")
                    print("[Protocols:Enzyme Solubility Error] There is something wrong with the [contact_number] section.")
                    quit()

                #Import our class
                try:
                    from pact.analysis.structure.contact_number import contact_number
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.structure.contact_number was not found.")
                    quit()
        
                #Create our object
                obj_contact = contact_number(self.obj_cfgparser, self.dict_programs, {})

                #Calculate the distance
                dict_contact = obj_contact.contact_number(dict_pdb)

            """
            *****************************************
            Output CSV and .pact
            *****************************************
            """
            quit()
            aa_table = 'ACDEFGHIKLMNPQRSTVWY'
            wtaa = self.obj_cfgparser.get('global', 'wtaa').upper()
            
            if self.dict_workflow['pdb_import']:
                chain = self.obj_cfgparser.get('classification_analysis', 'chain').upper()
                pdb_file = self.obj_cfgparser.get('classification_analysis', 'pdb_file')
                list_pdb_sites = sorted([x for x in dict_pdb[pdb_file]['dssp'][chain]])

            #Make a dict to work into
            dict_output = {}

            #Get the dataset name in order
            if self.obj_cfgparser.has_section("combinepact"):
                num_datasets = int(self.obj_cfgparser.get('combinepact', 'numdatasets'))
                list_datasets = [self.obj_cfgparser.get('combinepact', 'dataset_' + str(int_dataset)) 
                                 for int_dataset in range(1, num_datasets + 1)]
            else:
                list_datasets = []

            #Get the header
            str_output = ','.join([
                "Location",
                "Mutation",
                ','.join([dataset + "_fitness" for dataset in list_datasets]),
                ','.join([dataset + "_sd_from_wt" for dataset in list_datasets]),
                ','.join([dataset + "_classified" for dataset in list_datasets]),
                "wt_resi",
                "wt_pssm",
                "wt_percent",
                "max_pssm",
                "max_percent",
                "pssm_cons_count",
                "percent_cons_count",
                "wt_max_pssm",
                "wt_max_percent",
                "mut_pssm",
                "mut_percent",
                "frac_burial",
                "polarity",
                "aromatics",
                "philic_phobic",
                "size",
                "hydropathy",
                "dist_to_active",
                "contact_number"
                ]) + "\n"

            #Loop the locations
            for loc in range(1, len(wtaa) + 1):

                #Add to dict if not already added
                if loc not in dict_output:
                    dict_output[loc] = {}

                #Loop the mutations
                for mut in aa_table:

                    #Add to dict if not already added
                    if mut not in dict_output[loc]:
                        dict_output[loc][mut] = {}

                    #Get the location
                    str_output = str_output + str(loc) + ','

                    #Get the mutation
                    str_output = str_output + mut + ','

                    #Get the datasets, test if loc is in there
                    str_output = str_output + ','.join([str(dict_merged_datasets[dataset][loc][mut]['fitness'])
                                                        if loc in dict_merged_datasets[dataset] else " "
                                                        for dataset in list_datasets
                                                        ]) + ','

                    str_output = str_output + ','.join([str(dict_merged_datasets[dataset][loc][mut]['sd_from_wt'])
                                                        if loc in dict_merged_datasets[dataset] else " "
                                                        for dataset in list_datasets
                                                        ]) + ','

                    #Get the classified
                    str_output = str_output + ','.join([str(dict_classified[dataset][loc][mut])
                                                        if loc in dict_classified[dataset] else " "
                                                        for dataset in list_datasets
                                                        ]) + ','

                    #Get the WT Consensus Data
                    if self.dict_workflow['consensus']:
                        str_output = str_output + ','.join(map(str, [
                            dict_wtcons[loc]['wt_resi'],                                   
                            dict_wtcons[loc]['wt_pssm'],  
                            dict_wtcons[loc]['wt_percent'],  
                            dict_wtcons[loc]['max_pssm'],  
                            dict_wtcons[loc]['max_percent'],  
                            dict_wtcons[loc]['pssm_cons_count'],
                            dict_wtcons[loc]['percent_cons_count'],  
                            dict_wtcons[loc]['wt_max_pssm'],  
                            dict_wtcons[loc]['wt_max_percent'],  
                            ])) + ','
                    else:
                        str_output = str_output + "NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,"

                    #Get the PSSM data
                    if self.dict_workflow['pssm_reader']:
                        str_output = str_output + dict_pssm[loc][mut][0] + ","
                        str_output = str_output + dict_pssm[loc][mut][1] + ","
                    else:
                        str_output = str_output + "NaN,NaN,"

                    #Get the fraction burial
                    if self.dict_workflow['pdb_import']:
                        if loc in list_pdb_sites:
                            str_output = str_output + str(dict_pdb[pdb_file]['dssp'][chain][loc]['frac_burial']) + ","
                        else:
                            str_output = str_output + "NaN,"
                    else:
                        str_output = str_output + "NaN,"

                    #To/From Proline
                    if self.dict_workflow['residue_chemical_size']:
                        dict_rcs_mut = obj_rcs.mut_info(wtaa[loc - 1], mut)
                        str_output = str_output + ','.join([dict_rcs_mut['polarity'],
                                                            dict_rcs_mut['aromatics'],
                                                            dict_rcs_mut['philic_phobic'],
                                                            dict_rcs_mut['size'],
                                                            str(dict_rcs_mut['hydropathy']),
                                                            ]) + ','
                    else:
                        str_output = str_output + "NaN,NaN,NaN,NaN,NaN,"

                    #Dist to active site
                    if self.dict_workflow['distance_to_active']:
                        str_output = str_output + str(min(dict_dtoa_dist[chain][loc])) + ','
                    else:
                        str_output = str_output + "NaN,"
                    
                    #Contact number
                    if self.dict_workflow['contact_number']:
                        str_output = str_output + str(len(dict_contact[chain][loc])) + ','
                    else:
                        str_output = str_output + "NaN,"

                    #Newline
                    str_output = str_output + '\n'            

            #Output a csv file
            with open(self.directory + self.output_prefix + '_dataset.csv', 'w') as file_output:
                file_output.write(str_output)

            #At this point it's easier to backcalculate the csv file
            list_output = str_output.splitlines()
            list_keys = list_output[0].rstrip('\n').split(',')
            columns = len(list_keys)
            
            #Parse the lines
            for line in list_output[1:]:

                #Split the line
                splitline = line.split(',')

                #Parse each column
                for i in range(2, columns):
                    dict_output[int(splitline[0])][splitline[1]][list_keys[i]] = splitline[i]

            #Output a pact file
            print(save_pact_file(dict_output, self.directory + self.output_prefix + '_dataset'))

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
