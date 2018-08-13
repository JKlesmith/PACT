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

"""enzyme solubility - protocol to calculate mutations that pass our enzyme solubility filter"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Protocols:Enzyme Solubility Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from collections import Counter
from configparser import ConfigParser, NoSectionError, NoOptionError
from pact.pact_common import file_checker, pretty_counter_dicts
from sys import platform
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class enzyme_solubility:
    """Calculate mutations on a enzyme solubility filter"""

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
            print("[Protocols:Enzyme Solubility Error] The config file is incorrect.")
            print("[Protocols:Enzyme Solubility Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:Enzyme Solubility Error] The config file is incorrect.")
            print("[Protocols:Enzyme Solubility Error] There is something wrong with the [workflow] section (spelling???).")
            quit()

        #Get the PACT preamble
        self.pact_preamble = preamble

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Add trailing slash on directory
        if platform.startswith('win'):            
            if self.directory[-1:] != "\\":
                self.directory = self.directory + "\\"
        else:
            if self.directory[-1:] != "/":
                self.directory = self.directory + "/"
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
                print("[Protocols:Enzyme Solubility Error] The combinepact config file is incorrect.")
                print("[Protocols:Enzyme Solubility Error] There is something wrong with the [combinepact] section.")
                quit()

            #Import our combinepact class
            try:
                from pact.analysis.combine_pact import combine_pact
            except ImportError:
                print("[Protocols:Enzyme Solubility Error] combine_pact was not found.")

            #Create the object then call the merger
            obj_combine = combine_pact(self.obj_cfgparser, self.dict_programs, {})

            #The dict will be like {'dataset name': {data...
            dict_merged_datasets = obj_combine.combine_pact()

            """
            *****************************************
            Classify our mutations
            *****************************************
            """
            #Build a dict_classified with [location][mutation] = "DEL/SLIGHTDEL/NEU/NONE"
            #0.15 in GFP, 80% of WT = neutral, 50% of WT = slightly, <50% of WT = deleterious
            
            #Get the config file elements
            try:
                screen_dataset = self.obj_cfgparser.get("enzyme_solubility", "dataset_screen")
                screen_threshold = float(self.obj_cfgparser.get("enzyme_solubility", "screen_threshold"))
                fitness_dataset = self.obj_cfgparser.get("enzyme_solubility", "dataset_fitness")
                fitness_neutral = float(self.obj_cfgparser.get("enzyme_solubility", "fitness_neu"))
                fitness_slightdel = float(self.obj_cfgparser.get("enzyme_solubility", "fitness_slightdel"))
            except NoOptionError:
                print("[Enzyme Solubility Error] Missing [enzyme_solubility] config file elements.")
                quit()
            except ValueError:
                print("[Enzyme Solubility Error] Incorrect [enzyme_solubility] config file elements.")
                quit()
            except TypeError:
                print("[Enzyme Solubility Error] Incorrect [enzyme_solubility] config file elements.")
                quit()

            #Make a dict to add our classifications into
            dict_classified = {}
            dict_basal = {}

            #Loop the locations
            for loc in dict_merged_datasets[screen_dataset]:
                
                #Add a new location if not in the dict
                if loc not in dict_classified:
                    dict_classified[loc] = {}

                #Add a new location if not in the dict
                if loc not in dict_basal:
                    dict_basal[loc] = {}

                #Loop the muts
                for mut in dict_merged_datasets[screen_dataset][loc]:

                    #Skip WT, stop, and NaN
                    if (mut == dict_merged_datasets[fitness_dataset][loc][mut]['wt_residue'] or
                        mut == "*" or
                        dict_merged_datasets[fitness_dataset][loc][mut]['fitness'] == "NaN"):

                        dict_basal[loc][mut] = "UNCLASSIFIED"
                        continue

                    #Get the fitness value from the fitness dataset
                    fitness_value = float(dict_merged_datasets[fitness_dataset][loc][mut]['fitness'])

                    #For the basal screen fitness
                    #Assign a classification of deleterious, slightly deleterious, or neutral
                    if fitness_value < fitness_slightdel:
                        dict_basal[loc][mut] = "DEL"

                    elif (fitness_value >= fitness_slightdel and fitness_value < fitness_neutral):
                        dict_basal[loc][mut] = "SLIGHTDEL"

                    elif fitness_value >= fitness_neutral:
                        dict_basal[loc][mut] = "NEU"

                    #Skip WT, stop, and NaN
                    if (mut == dict_merged_datasets[screen_dataset][loc][mut]['wt_residue'] or
                        mut == "*" or
                        dict_merged_datasets[screen_dataset][loc][mut]['fitness'] == "NaN"):

                        dict_classified[loc][mut] = "UNCLASSIFIED"
                        continue

                    #Are we are enriched in the screen dataset?
                    if float(dict_merged_datasets[screen_dataset][loc][mut]['fitness']) < screen_threshold:
                        dict_classified[loc][mut] = "UNCLASSIFIED"
                        continue

                    #Assign a classification of deleterious, slightly deleterious, or neutral
                    if fitness_value < fitness_slightdel:
                        dict_classified[loc][mut] = "DEL"

                    elif (fitness_value >= fitness_slightdel and fitness_value < fitness_neutral):
                        dict_classified[loc][mut] = "SLIGHTDEL"

                    elif fitness_value >= fitness_neutral:
                        dict_classified[loc][mut] = "NEU"

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
                    print("[Protocols:Enzyme Solubility] pact.analysis.basal_count was not found.")
                    quit()
        
                #Create our object
                obj_basal = basal_count(self.obj_cfgparser, self.dict_programs, {})

                #Count our basal rates
                print("[Protocols:Enzyme Solubility] Basal Screen Counts")
                file_output.write("[Protocols:Enzyme Solubility] Basal Screen Counts")
                file_output.write(obj_basal.basal_count(dict_basal) + "\n")

                print("[Protocols:Enzyme Solubility] Basal Fitness Counts")
                file_output.write("[Protocols:Enzyme Solubility] Basal Fitness Counts")
                file_output.write(obj_basal.basal_count(dict_classified) + "\n")

            """
            *****************************************
            DNA Filtering/Alignment or PSSM Object
            *****************************************
            """
            if self.dict_workflow['blastp_align_filter'] or self.dict_workflow['pssm']:

                #Import Check Output
                from subprocess import check_output

                #Import our class
                try:
                    from pact.analysis.sequence.homology_pssm import homology_classifier
                except ImportError:
                    print("[Protocols:Enzyme Solubility] pact.analysis.sequence.homology_pssm was not found.")
                    quit()
        
                #Create our object
                obj_homology = homology_classifier(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

            """
            *****************************************
            DNA Filtering/Alignment
            *****************************************
            """
            if self.dict_workflow['blastp_align_filter']:
               
                #Convert our XML file
                print("[Protocols:Enzyme Solubility] xml_to_fasta")
                file_output.write("[Protocols:Enzyme Solubility] xml_to_fasta\n")
                obj_homology.xml_to_fasta()

                #Run CD-HIT on our new fasta file
                print("[Protocols:Enzyme Solubility] cdhit")
                file_output.write("[Protocols:Enzyme Solubility] cdhit\n")

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
                print("[Protocols:Enzyme Solubility] cdhit_wtcheck")
                file_output.write("[Protocols:Enzyme Solubility] cdhit_wtcheck\n")
                obj_homology.cdhit_wt_check()

                #Run MUSCLE on our new fasta file
                print("[Protocols:Enzyme Solubility] muscle")
                file_output.write("[Protocols:Enzyme Solubility] muscle\n")
                check_output([self.dict_programs['muscle'],
                              "-in",
                              self.directory + self.output_prefix + ".afa",
                              "-out",
                              self.directory + self.output_prefix + ".msa"])

                #Process our MSA (needs to be on for PSIBlast)
                print("[Protocols:Enzyme Solubility] processmsa")
                file_output.write("[Protocols:Enzyme Solubility] processmsa\n")
                list_msa = obj_homology.process_msa()

                #Save our list
                print("[Protocols:Enzyme Solubility] Saving our MSA")
                file_output.write("[Protocols:Enzyme Solubility] Saving our MSA\n")
                obj_homology.save_data_structure(list_msa, "list_msa")

            """
            *****************************************
            PSSM
            *****************************************
            """
            if self.dict_workflow['pssm'] or self.dict_workflow['strict_filter']:
                #Open our list
                print("[Protocols:Enzyme Solubility] Opening our MSA")
                file_output.write("[Protocols:Enzyme Solubility] Opening our MSA\n")
                list_msa = obj_homology.open_data_structure("list_msa")

                #Split our msa for PSIBlast (needs to be on for PSIBlast)
                print("[Protocols:Enzyme Solubility] msa_split")
                file_output.write("[Protocols:Enzyme Solubility] msa_split\n")
                list_pbcmds = obj_homology.msa_split(list_msa)

                #Run PSIBlast
                print("[Protocols:Enzyme Solubility] psiblast")
                file_output.write("[Protocols:Enzyme Solubility] psiblast\n")
                for command in list_pbcmds:
                    check_output([self.dict_programs['psiblast'], *command])

                #Import our PSSM data
                print("[Protocols:Enzyme Solubility] pssm_file_import")
                file_output.write("[Protocols:Enzyme Solubility] pssm_file_import\n")
                dict_pssm = obj_homology.pssm_file_import()

                #Save our heatmap
                print("[Protocols:Enzyme Solubility] Saving a PSSM .csv heatmap")
                file_output.write(obj_homology.pssm_output_heat(dict_pssm) + "\n")

                #Save our csv
                print("[Protocols:Enzyme Solubility] Saving a PSSM .csv column data")
                file_output.write(obj_homology.pssm_output_csv(dict_pssm) + "\n")

                #Save our PACT File
                print("[Protocols:Enzyme Solubility] Saving a PSSM .pact file")
                file_output.write(obj_homology.save_data_structure(dict_pssm, "PSSM") + "\n")

                #Count our classifiers
                print("[Protocols:Enzyme Solubility] PSSM")
                file_output.write("[Protocols:Enzyme Solubility] PSSM")

                print("Fitness Rates:")
                file_output.write("Fitness Rates:")
                file_output.write(obj_homology.classified_count_pssm(dict_pssm, dict_basal) + "\n")
                
                print("\nScreen Rates:")
                file_output.write("\nScreen Rates:")
                file_output.write(obj_homology.classified_count_pssm(dict_pssm, dict_classified) + "\n")

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

                #Return the process dict {1: {'A': {''
                dict_rcs = obj_rcs.process_dataset(dict_merged_datasets)

                #Count our classifiers
                print("[Protocols:Enzyme Solubility] Residue Chemical/Size")
                file_output.write("[Protocols:Enzyme Solubility] Residue Chemical/Size")

                print("Fitness Rates:")
                file_output.write("Fitness Rates:")
                file_output.write(obj_rcs.classified_count(dict_rcs, fitness_dataset, dict_basal) + "\n")
                
                print("\nScreen Rates:")
                file_output.write("\nScreen Rates:")
                file_output.write(obj_rcs.classified_count(dict_rcs, screen_dataset, dict_classified) + "\n")

            """
            *****************************************
            PDB Import Section
            *****************************************
            """
            #Only import and run if selected
            if (self.dict_workflow['pdb_import'] or 
                self.dict_workflow['distance_to_active'] or 
                self.dict_workflow['contact_number'] or 
                self.dict_workflow['strict_filter']):

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('pdb_import'):           
                    print("[Protocols:Enzyme Solubility Error] The pdb_import config file is incorrect.")
                    print("[Protocols:Enzyme Solubility Error] There is something wrong with the [pdb_import] section.")
                    quit()

                #Import our combinepact class
                try:
                    from pact.analysis.pdb_import import pdb_import
                except ImportError:
                    print("[Protocols:Enzyme Solubility Error] pdb_import was not found.")

                #Create the object then call the merger
                obj_pdb = pdb_import(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #The dict will be like {'pdb name': {data...
                dict_pdb = obj_pdb.pdb_import()

            """
            *****************************************
            Distance to Active Site
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['distance_to_active'] or self.dict_workflow['strict_filter']:

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

                #Count our classifiers
                print("[Protocols:Enzyme Solubility] Distance to Active Site")
                file_output.write("[Protocols:Enzyme Solubility] Distance to Active Site")

                print("Fitness Rates:")
                file_output.write("Fitness Rates:")
                file_output.write(obj_dtoa.classified_count(dict_dtoa_dist, dict_basal) + "\n")
                
                print("\nScreen Rates:")
                file_output.write("\nScreen Rates:")
                file_output.write(obj_dtoa.classified_count(dict_dtoa_dist, dict_classified) + "\n")

            """
            *****************************************
            Contact Number
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['contact_number'] or self.dict_workflow['strict_filter']:

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

                #Count our classifiers
                print("[Protocols:Enzyme Solubility] Contact Number")
                file_output.write("[Protocols:Enzyme Solubility] Contact Number")

                print("Fitness Rates:")
                file_output.write("Fitness Rates:")
                file_output.write(obj_contact.classified_count(dict_contact, dict_basal) + "\n")
                
                print("\nScreen Rates:")
                file_output.write("\nScreen Rates:")
                file_output.write(obj_contact.classified_count(dict_contact, dict_classified) + "\n")

            """
            *****************************************
            Strict Enzyme Filter
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['strict_filter']:

                print("[Protocols:Enzyme Solubility] Strict Enzyme Filter")
                file_output.write("[Protocols:Enzyme Solubility] Strict Enzyme Filter")

                #Check if the dicts exist
                if 'dict_pssm' not in locals():
                    print("[Protocols:Enzyme Solubility] Missing PSSM Data")
                    file_output.write("[Protocols:Enzyme Solubility] Missing PSSM Data")
                    quit()

                if 'dict_contact' not in locals():
                    print("[Protocols:Enzyme Solubility] Missing Contact Number Data")
                    file_output.write("[Protocols:Enzyme Solubility] Missing Contact Number Data")
                    quit()

                if 'dict_dtoa_dist' not in locals():
                    print("[Protocols:Enzyme Solubility] Missing Active Site Distance Data")
                    file_output.write("[Protocols:Enzyme Solubility] Missing Active Site Distance Data")
                    quit()


                #Implement a filter that
                #PSSM >= 0
                #Distance to Active >= 15A
                #Contact Number <= 16
                #No proline mutations

                #Create a list to work into
                list_strictfilter = []

                #Loop the locations
                for loc in dict_classified:

                    #Loop the mutations
                    for mut in dict_classified[loc]:

                        #Skip PRO, and Stop
                        if (mut == "P" or 
                            mut == "*" or 
                            dict_merged_datasets[screen_dataset][loc][mut]['wt_residue'] == "P"):
                            continue

                        #Check if PSSM is less than 0
                        if int(dict_pssm[loc][mut][0]) < 0:
                            continue

                        #Skip residues without location data
                        if loc not in dict_contact or loc not in dict_dtoa_dist:
                            continue

                        #Check if Distance to Active is less than 15A
                        if min(dict_dtoa_dist[loc]) < 15:
                            continue

                        #Check if the contact number is greater than 16
                        if len(dict_contact[loc]) > 16:
                            continue

                        #Otherwise, add to our list
                        list_strictfilter.append(dict_classified[loc][mut])

                #Report
                str_return = '\n'.join(map(str, [
                "Enzyme Strict Filter",
                pretty_counter_dicts(dict(Counter(list_strictfilter)))
                ]))

                print(str_return)
                file_output.write(str_return)

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:Enzyme Solubility Error] This protocol needs to be ran within the context of PACT.")
