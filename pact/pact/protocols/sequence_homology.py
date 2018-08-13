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

"""sequence homology - compare sequence homology to our selections"""

from sys import version_info

#Setup our protocol name
str_protocol_name = "Sequence Homology"

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
from subprocess import check_output
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class sequence_homology:
    """Calculate the sequence homology and PSSM"""

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
            print("[Protocols:Homology Error] The config file is incorrect.")
            print("[Protocols:Homology Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:Homology Error] The config file is incorrect.")
            print("[Protocols:Homology Error] There is something wrong with the [workflow] section (spelling???).")
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

            #Import our class
            try:
                from pact.analysis.sequence.homology_pssm import homology_classifier
            except ImportError:
                print("[Protocols:Homology Error] pact.analysis.sequence.homology_pssm was not found.")
                quit()
        
            #Create our object
            obj_homology = homology_classifier(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

            """
            *****************************************
            DNA Filtering/Alignment Section
            *****************************************
            """
            if self.dict_workflow['blastp_align_filter']:
               
                #Convert our XML file
                print("[Protocols:Homology] xml_to_fasta")
                file_output.write("[Protocols:Homology] xml_to_fasta\n")
                obj_homology.xml_to_fasta()

                #Run CD-HIT on our new fasta file
                print("[Protocols:Homology] cdhit")
                file_output.write("[Protocols:Homology] cdhit\n")

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
                print("[Protocols:Homology] cdhit_wtcheck")
                file_output.write("[Protocols:Homology] cdhit_wtcheck\n")
                obj_homology.cdhit_wt_check()

                #Run MUSCLE on our new fasta file
                print("[Protocols:Homology] muscle")
                file_output.write("[Protocols:Homology] muscle\n")
                check_output([self.dict_programs['muscle'],
                              "-in",
                              self.directory + self.output_prefix + ".afa",
                              "-out",
                              self.directory + self.output_prefix + ".msa"])

                #Process our MSA (needs to be on for PSIBlast)
                print("[Protocols:Homology] processmsa")
                file_output.write("[Protocols:Homology] processmsa\n")
                list_msa = obj_homology.process_msa()

                #Save our list
                print("[Protocols:Homology] Saving our MSA")
                file_output.write("[Protocols:Homology] Saving our MSA\n")
                save_pact_file(list_msa, self.directory + self.output_prefix + '_' + "list_msa")

            """
            *****************************************
            PSSM Section
            *****************************************
            """
            if self.dict_workflow['pssm']:
                #Open our list
                print("[Protocols:Homology] Opening our MSA")
                file_output.write("[Protocols:Homology] Opening our MSA\n")
                list_msa = open_pact_file(self.directory + self.output_prefix + '_' + "list_msa")

                #Split our msa for PSIBlast (needs to be on for PSIBlast)
                print("[Protocols:Homology] msa_split")
                file_output.write("[Protocols:Homology] msa_split\n")
                list_pbcmds = obj_homology.msa_split(list_msa)

                #Run PSIBlast
                print("[Protocols:Homology] psiblast")
                file_output.write("[Protocols:Homology] psiblast\n")
                for command in list_pbcmds:
                    check_output([self.dict_programs['psiblast'], *command])

                #Import our PSSM data
                print("[Protocols:Homology] pssm_file_import")
                file_output.write("[Protocols:Homology] pssm_file_import\n")
                dict_pssm = obj_homology.pssm_file_import()

                #Save our heatmap
                print("[Protocols:Homology] Saving a PSSM .csv heatmap")
                file_output.write(obj_homology.pssm_output_heat(dict_pssm) + "\n")

                #Save our csv
                print("[Protocols:Homology] Saving a PSSM .csv column data")
                file_output.write(obj_homology.pssm_output_csv(dict_pssm) + "\n")

                #Save our PACT File
                print("[Protocols:Homology] Saving a PSSM .pact file")
                file_output.write(save_pact_file(dict_pssm, self.directory + self.output_prefix + '_' + "PSSM") + "\n")

            """
            *****************************************
            Sitewise Frequency Section
            *****************************************
            """
            if self.dict_workflow['site_frequencies']:
                #Open our list
                print("[Protocols:Homology] Opening our MSA")
                file_output.write("[Protocols:Homology] Opening our MSA\n")
                list_msa = open_pact_file(self.directory + self.output_prefix + '_' + "list_msa")

                #Calculate our frequencies
                print("[Protocols:Homology] Calculate our frequencies")
                file_output.write("[Protocols:Homology] Calculate our frequencies\n")
                dict_freq = obj_homology.msa_freq(list_msa)

                #Save our CSV heatmap
                print("[Protocols:Homology] Saving the frequencies heatmap")
                file_output.write("[Protocols:Homology] Saving the frequencies heatmap\n")
                obj_homology.freq_output_heat(dict_freq)

                #Save our PACT File
                print("[Protocols:Homology] Saving a Freq .pact file")
                file_output.write(save_pact_file(dict_freq, self.directory + self.output_prefix + '_' + "freq") + "\n")

            """
            *****************************************
            Read stored .pact files
            *****************************************
            """
            if self.dict_workflow['pssm_reader']:
                #Open our PACT File
                print("[Protocols:" + str_protocol_name + "] Opening a PSSM .pact file")
                dict_pssm = open_pact_file(self.directory + self.output_prefix + '_' + "PSSM")

                #Count our classifiers
                print("[Protocols:" + str_protocol_name + "] PSSM Classifier Count")
                file_output.write("[Protocols:" + str_protocol_name + "] PSSM Classifier Count")

            """
            *****************************************
            Pact Combine Section
            *****************************************
            """
            if self.dict_workflow['combinepact']:

                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('combinepact'):           
                    print("[Protocols:Homology Error] The combinepact config file is incorrect.")
                    print("[Protocols:Homology Error] There is something wrong with the [combinepact] section.")
                    quit()

                #Import our combinepact class
                try:
                    from pact.analysis.combine_pact import combine_pact
                except ImportError:
                    print("[Protocols:Homology Error] combine_pact was not found.")

                #Create the object then call the merger
                obj_combine = combine_pact(self.obj_cfgparser, self.dict_programs, {})

                #The dict will be like {'dataset name': {data...
                dict_merged_datasets = obj_combine.combine_pact()

            """
            *****************************************
            Analysis Section
            *****************************************
            """
            if self.dict_workflow['analysis_sitefitness_homology']:
                #Which dataset do we want?
                if self.obj_cfgparser.get('analysis_sitefitness_homology', 'dataset_x') == "site_frequencies":
                    file_name = "freq"                    
                else:
                    file_name = "pssm"

                #Open our dict
                print("[Protocols:Homology] Opening our homology data")
                file_output.write("[Protocols:Homology] Opening our site freqs\n")
                dict_homology = open_pact_file(self.directory + self.output_prefix + '_' + file_name)

                #Plot our data
                if self.obj_cfgparser.get('analysis_sitefitness_homology', 'scatter') == "True":
                    print("[Protocols:Homology] Plotting the figure")
                    file_output.write("[Protocols:Homology] Plotting the figure\n")
                    obj_homology.analysis_site_fit_homology_plot(dict_homology, dict_merged_datasets, file_name)

                #Plot our data
                print("[Protocols:Homology] Making our classifier table")
                file_output.write("[Protocols:Homology] Making our classifier table\n")
                obj_homology.analysis_site_fit_homology_classifier(dict_homology, dict_merged_datasets, file_name)

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
