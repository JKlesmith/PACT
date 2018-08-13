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

"""back to consensus - For a given alignment of homologous sequences, what is the hit rate from the selection or shared selections"""

from sys import version_info

#Setup our protocol name
str_protocol_name = "Back to Consensus"

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Protocols:" + str_protocol_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
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

class back_to_consensus:
    """For a given alignment of homologous sequences, what is the hit rate from the selection or shared selections"""

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
            Pact Combine (Required)
            *****************************************
            """
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
            if self.dict_workflow['pssm_reader'] or self.dict_workflow['wt_consensus']:
                #Open our PACT File
                print("[Protocols:" + str_protocol_name + "] Opening a PSSM .pact file")
                dict_pssm = open_pact_file(self.directory + self.output_prefix + '_' + "PSSM")

                #Count our classifiers
                print("[Protocols:" + str_protocol_name + "] PSSM Classifier Count")
                file_output.write("[Protocols:" + str_protocol_name + "] PSSM Classifier Count")
                

                for dataset in dict_classified:
                    print("[Protocols:" + str_protocol_name + "] PSSM Fitness Rates for dataset: " + dataset)
                    file_output.write("[Protocols:" + str_protocol_name + "] PSSM Fitness Rates for dataset: " + dataset)
                    file_output.write(obj_homology.classified_count_pssm(dict_pssm, dict_classified[dataset]) + "\n")

                    print("[Protocols:" + str_protocol_name + "] Wrote CSV of fitness values categorized by PSSM group and mutation type for dataset: " + dataset)
                    file_output.write("[Protocols:" + str_protocol_name + "] PSSM Fitness Rates for dataset: " + dataset)
                    file_output.write(obj_homology.classified_count_pssm_csv(dict_pssm, dict_classified[dataset], 
                                                                         dict_merged_datasets, dataset, class_column) + "\n")

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

                #Import our combinepact class
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

                #Get the prob of finding a classified mutation
                for dataset in dict_classified:
                    print("[Protocols:" + str_protocol_name + "] Mut Classification vs WT Cons for dataset: " + dataset)
                    obj_consensus.wtcons_count_class(dict_wtcons, dict_classified[dataset], dataset)

                #Get the prob of finding a classified mutation
                for dataset in dict_classified:
                    print("[Protocols:" + str_protocol_name + "] Mutating a non-conserved site to a conserved site: " + dataset)
                    file_output.write(obj_consensus.nonconserved_sites(dict_wtcons, dict_pssm, dict_classified[dataset], dataset))

                #Get the prob of finding a classified mutation
                for dataset in dict_classified:
                    print("[Protocols:" + str_protocol_name + "] Mutating a non-conserved site to a non-conserved mutation: " + dataset)
                    file_output.write(obj_consensus.nonconserved_mutations(dict_wtcons, dict_pssm, dict_classified[dataset], dataset))

                #Get the cross set distribution
                obj_consensus.cons_count_setvset(dict_wtcons, dict_pssm, dict_classified)

                if self.dict_workflow['pdb_import']:
                    for dataset in dict_classified:
                        print("[Protocols:" + str_protocol_name + 
                              "] Mutating a non-conserved site to a conserved site (Buried Residues Only): " + dataset)
                        file_output.write(obj_consensus.nonconserved_sites_burial(
                            dict_wtcons, dict_pssm, dict_classified[dataset], dataset, dict_pdb, "<"))

                    for dataset in dict_classified:
                        print("[Protocols:" + str_protocol_name + 
                              "] Mutating a non-conserved site to a conserved site (Surface Residues Only): " + dataset)
                        file_output.write(obj_consensus.nonconserved_sites_burial(
                            dict_wtcons, dict_pssm, dict_classified[dataset], dataset, dict_pdb, ">="))

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
