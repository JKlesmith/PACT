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

"""tools - call our tool scripts"""

from sys import version_info

#Setup our protocol name
str_protocol_name = "Tools"

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
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class tools:
    """Provide an entrypoint for the tools scripts"""

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
            print("[Protocols:Tools Error] The config file is incorrect.")
            print("[Protocols:Tools Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:Tools Error] The config file is incorrect.")
            print("[Protocols:Tools Error] There is something wrong with the [workflow] section (spelling???).")
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
            fastq_split Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['fastq_split']:

                #Set custom locations for the fastq files
                if self.obj_cfgparser.get('fastq_split', 'directory') != "":
                    fastq_dir = self.obj_cfgparser.get('fastq_split', 'directory')
                else:
                    fastq_dir = self.directory

                #Do a check if the user accidently lists the same file for the fwd and rev
                if self.obj_cfgparser.get('fastq_split', 'forward_fastq') == self.obj_cfgparser.get('fastq_split', 'reverse_fastq'):
                    print("[Protocols:Tools Error] The selected forward and reverse fastq files are the same.")
                    quit()

                #Create our options dicts
                try:
                    dict_fastq_split = {
                        'forward_file':fastq_dir + self.obj_cfgparser.get('fastq_split', 'forward_fastq'),
                        'reverse_file':fastq_dir + self.obj_cfgparser.get('fastq_split', 'reverse_fastq'),
                        'output_prefix':self.directory,
                        'cutoff':self.obj_cfgparser.get('fastq_split', 'cutoff')
                        }
                except NoSectionError:
                    print("[Protocols:Tools Error] The fastq_split config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the [fastq_split_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Tools Error] The fastq_split config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the name of a option flag.")
                    quit()
            
                #Import our fastq_split file
                try:
                    from pact.tools.fastq_split import fastq_split
            
                except ImportError:
                    print("[Protocols:Tools Error] fastq_split was not found.")

                #Create the object
                obj_fastq_split = fastq_split(dict_fastq_split)
        
                #Call the entrypoint
                print("[Protocols:Tools] Splitting the fastq files.")
                file_output.write("[Protocols:Tools] Splitting the fastq files.\n")
                file_output.write(obj_fastq_split.fastq_split() + "\n")

            """
            *****************************************
            codon_condenser Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['codon_condenser']:
            
                #Create our options dicts
                try:
                    list_aminos = self.obj_cfgparser.get('codon_condenser', 'list_aminos')
                    codon = self.obj_cfgparser.get('codon_condenser', 'codon')
                except NoSectionError:
                    print("[Protocols:Tools Error] The codon_condenser config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the [codon_condenser_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Tools Error] The codon_condenser config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the name of a option flag.")
                    quit()
            
                #Import our fastq_split file
                try:
                    from pact.tools.codon_condenser import codon_condenser
            
                except ImportError:
                    print("[Protocols:Tools Error] codon_condenser was not found.")

                #Create the object then call the first file
                obj_cc = codon_condenser()
        
                if len(codon) == 3:
                    print("[Protocols:Tools] Codon analysis.")
                    file_output.write("[Protocols:Tools] Codon analysis.\n")
                    dict_analysis, ca_output = obj_cc.codon_analysis(codon)
                    print(ca_output)
                    file_output.write(ca_output + "\n")

                if len(list_aminos) > 0:
                    print("[Protocols:Tools] Codon condense.")
                    file_output.write("[Protocols:Tools] Codon condense.\n")
                    dict_analysis, cc_output = obj_cc.codon_condense(list_aminos)
                    print(cc_output)
                    file_output.write(cc_output + "\n")

            """
            *****************************************
            codon_swap Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['codon_swap']:
                #Create our options dicts
                try:
                    dict_options = {
                    'WT_DNA':self.obj_cfgparser.get('codon_swap', 'dna_sequence')
                    }
                except NoSectionError:
                    print("[Protocols:Tools Error] The codon_swap config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the [codon_swap_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Tools Error] The codon_swap config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the name of a option flag.")
                    quit()
            
                #Import our fastq_split file
                try:
                    from pact.tools.codon_swap import codon_swap
            
                except ImportError:
                    print("[Protocols:Tools Error] codon_condenser was not found.")

                #Create the object then call the first file
                obj_cs = codon_swap(dict_options)

                print("[Protocols:Tools] Codon swap.")
                file_output.write("[Protocols:Tools] Codon swap.\n")
                new_dna_seq, cs_output = obj_cs.codon_swap()
                print(cs_output)
                file_output.write(cs_output + "\n")

            """
            *****************************************
            fastq_to_fasta Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['fastq_to_fasta']:
                #Create our options dicts
                try:
                    dict_options = {
                    'fastq_file':self.directory + self.obj_cfgparser.get('fastq_to_fasta', 'fastq_file')
                    }
                except NoSectionError:
                    print("[Protocols:Tools Error] The fastq_to_fasta config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the [fastq_to_fasta_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Tools Error] The fastq_to_fasta config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the name of a option flag.")
                    quit()
            
                #Import our fastq_split file
                try:
                    from pact.tools.fastq_to_fasta import fastq_to_fasta
            
                except ImportError:
                    print("[Protocols:Tools Error] fastq_to_fasta was not found.")

                #Create the object then call the first file
                obj_ftf = fastq_to_fasta(dict_options)

                print("[Protocols:Tools] fastq_to_fasta.")
                file_output.write("[Protocols:Tools] fastq_to_fasta.\n")
                obj_ftf.fastq_to_fasta()

            """
            *****************************************
            primer_design Section
            *****************************************
            """
            #Only import and run if selected
            if self.dict_workflow['primer_design']:
                #Create our options dicts
                try:
                    dict_options = {
                    'WTDNA':self.obj_cfgparser.get('primer_design', 'dna_sequence'),
                    'Processes':self.obj_cfgparser.get('primer_design', 'processes'),
                    'Out_Prefix':self.directory + self.output_prefix,
                    'mutcodons':self.obj_cfgparser.get('primer_design', 'mutated_codons')
                    }
                except NoSectionError:
                    print("[Protocols:Tools Error] The primer_design config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the [primer_design_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Tools Error] The primer_design config file is incorrect.")
                    print("[Protocols:Tools Error] There is something wrong with the name of a option flag.")
                    quit()
            
                #Import our fastq_split file
                try:
                    from pact.tools.primer_design import primer_design
            
                except ImportError:
                    print("[Protocols:Tools Error] codon_condenser was not found.")

                #Create the object then call the first file
                obj_pd = primer_design(dict_options)

                print("[Protocols:Tools] Primer design.")
                file_output.write("[Protocols:Tools] Primer design.\n")
                obj_pd.primer_design()

            """
            *****************************************
            Read in fitness csv
            *****************************************
            """
            if self.dict_workflow['convert_csv_to_pact']:
                #Check to see if the section is there
                if not self.obj_cfgparser.has_section('convert_csv_to_pact'):           
                    print("[Protocols:" + str_protocol_name + " Error] The convert_csv_to_pact config file is incorrect.")
                    print("[Protocols:" + str_protocol_name + " Error] There is something wrong with the [convert_csv_to_pact] section.")
                    quit()

                #Import our convert_csv_to_pact class
                try:
                    from pact.analysis.convert_csv_to_pact import convert_csv_to_pact
                except ImportError:
                    print("[Protocols:" + str_protocol_name + " Error] convert_csv_to_pact was not found.")

                #Create the object then call the merger
                obj_csv = convert_csv_to_pact(self.obj_cfgparser, self.dict_programs, {'directory':self.directory})

                #The dict will be like {'dataset name': {data...
                obj_csv.read_csv_fitness()

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
