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

"""fitness - pipeline run to go from raw fastq reads to a fitness metric"""

from sys import version_info

#Setup our protocol name and version
str_protocol_name = "Fitness"

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Protocols:" + str_protocol_name + " Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from ast import literal_eval
from os import access, W_OK, R_OK
from configparser import ConfigParser, NoSectionError, NoOptionError
from pact.pact_common import file_checker
from time import strftime

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class fitness:

    def __init__(self, file_config, dict_programs, preamble):
        """Assign the class variables from our config file."""

        #Read our config file
        obj_cfgparser = ConfigParser()
        obj_cfgparser.read(file_config)

        #Send the cfgparser object as a global class object
        self.obj_cfgparser = obj_cfgparser
        
        #Assign our class variables from the config file
        self.wtdna = obj_cfgparser.get("global", "wtdna").upper()
        self.wtaa = obj_cfgparser.get("global", "wtaa").upper()
        self.mutcodons = obj_cfgparser.get('global', 'mutcodons')
        self.mutationtype = obj_cfgparser.get('global', 'mutationtype').lower()
        self.output_prefix = obj_cfgparser.get("global", "output_prefix")
        self.directory = obj_cfgparser.get("global", "directory")

        #Test the mutation type
        if self.mutationtype != "single" and self.mutationtype != "multiple":
            self.mutationtype = "single"

        #Assign these but catch value errors
        try:
            self.firstaamutated = int(obj_cfgparser.get("global", "firstaamutated"))
            self.lastaamutated = int(obj_cfgparser.get("global", "lastaamutated"))
            self.processes = int(obj_cfgparser.get("global", "processes"))
            self.mutthreshold = int(obj_cfgparser.get('global', 'mutthreshold'))
        except ValueError:
            print("Incorrect value for the global settings (string entered when it should be a number).")
            quit()

        #Read in the workflow
        try:
            self.dict_workflow = {mapping[0].lower(): literal_eval(mapping[1]) 
                              for mapping in obj_cfgparser.items("workflow")}
        except NoSectionError:
            print("[Protocols:Fitness Error] The config file is incorrect.")
            print("[Protocols:Fitness Error] There is something wrong with the [workflow] section.")
            quit()
        except ValueError:
            print("[Protocols:Fitness Error] The config file is incorrect.")
            print("[Protocols:Fitness Error] There is something wrong with the [workflow] section (spelling???).")
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
            FASTQ Merge Section
            *****************************************
            """
            if self.dict_workflow['fastq_merge_sel'] or self.dict_workflow['fastq_merge_ref']:
                #Import our fastq_merge file
                try:
                    from pact.sequencing.fastq_merge import fastq_merge
                except ImportError:
                    print("[Protocols:Fitness Error] fastq_merge was not found.")

            if self.dict_workflow['fastq_merge_ref']:
                #Set custom locations for the fastq files
                if self.obj_cfgparser.get('fastq_merge_ref', 'directory') != "":
                    ref_fastq_dir = self.obj_cfgparser.get('fastq_merge_ref', 'directory')
                else:
                    ref_fastq_dir = self.directory

                if self.obj_cfgparser.get('fastq_merge_ref', 'forward_fastq') == self.obj_cfgparser.get('fastq_merge_ref', 'reverse_fastq'):
                    print("[Protocols:Fitness Error] The reference forward and reverse fastq files are the same.")
                    quit()

                #Create our options dicts
                try:
                    dict_fastqmerge_options_ref = {
                        'Processes':self.processes,
                        'Forward_FASTQ':ref_fastq_dir + self.obj_cfgparser.get('fastq_merge_ref', 'forward_fastq'),
                        'Reverse_FASTQ':ref_fastq_dir + self.obj_cfgparser.get('fastq_merge_ref', 'reverse_fastq'),
                        'Out_Prefix':self.directory + self.output_prefix + "_Ref",
                        'Min_Coverage':float(self.obj_cfgparser.get('fastq_merge_ref', 'min_coverage')),
                        }
                except NoSectionError:
                    print("[Protocols:Fitness Error] The fastq_merge config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [fastq_merge_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The fastq_merge config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()

                #Create the object then call the first file
                obj_fastqmerge_ref = fastq_merge(dict_fastqmerge_options_ref)
       
                print("[Protocols:Fitness] Merging the reference fastq files.")
                file_output.write("[Protocols:Fitness] Merging the reference fastq files.\n")
                file_output.write(obj_fastqmerge_ref.fastq_merge() + "\n")

            if self.dict_workflow['fastq_merge_sel']:
                #Set custom locations for the fastq files
                if self.obj_cfgparser.get('fastq_merge_sel', 'directory') != "":
                    sel_fastq_dir = self.obj_cfgparser.get('fastq_merge_sel', 'directory')
                else:
                    sel_fastq_dir = self.directory

                #Do a check if the user accidently lists the same file for the fwd and rev
                if self.obj_cfgparser.get('fastq_merge_sel', 'forward_fastq') == self.obj_cfgparser.get('fastq_merge_sel', 'reverse_fastq'):
                    print("[Protocols:Fitness Error] The selected forward and reverse fastq files are the same.")
                    quit()

                #Create our options dicts
                try:
                    dict_fastqmerge_options_sel = {
                        'Processes':self.processes,
                        'Forward_FASTQ':sel_fastq_dir + self.obj_cfgparser.get('fastq_merge_sel', 'forward_fastq'),
                        'Reverse_FASTQ':sel_fastq_dir + self.obj_cfgparser.get('fastq_merge_sel', 'reverse_fastq'),
                        'Out_Prefix':self.directory + self.output_prefix + "_Sel",
                        'Min_Coverage':float(self.obj_cfgparser.get('fastq_merge_sel', 'min_coverage')),
                        }
                except NoSectionError:
                    print("[Protocols:Fitness Error] The fastq_merge config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [fastq_merge_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The fastq_merge config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()
            
                #Create the object then call the first file
                obj_fastqmerge_sel = fastq_merge(dict_fastqmerge_options_sel)
        
                print("[Protocols:Fitness] Merging the selected fastq files.")
                file_output.write("[Protocols:Fitness] Merging the selected fastq files.\n")
                file_output.write(obj_fastqmerge_sel.fastq_merge() + "\n")

            """
            *****************************************
            FASTQ Filter Section
            *****************************************
            """
            if self.dict_workflow['fastq_filter_translate_ref'] or self.dict_workflow['fastq_filter_translate_sel']:
                #Import our fastq_reader file
                try:
                    from pact.sequencing.fastq_filter_translate import fastq_filter_translate
                except ImportError:
                    print("[Protocols:Fitness Error] fastq_filter_translate was not found.")

            if self.dict_workflow['fastq_filter_translate_ref']:
            
                try:
                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('fastq_filter_translate_ref', 'fastq_file')) > 0:
                        file_input_fastq_ref = self.obj_cfgparser.get('fastq_filter_translate_ref', 'fastq_file')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Ref_Merge.fastq"):
                            file_input_fastq_ref = self.directory + self.output_prefix + "_Ref_Merge.fastq"

                    dict_fastqread_options_ref = {
                    'WTDNA':self.wtdna,
                    'WTAA':self.wtaa,
                    'FirstAAMutated':self.firstaamutated,
                    'LastAAMutated':self.lastaamutated,
                    'Processes':self.processes,
                    '5pAnchor':self.obj_cfgparser.get('fastq_filter_translate_ref', 'fiveprimeanchor'),
                    'MutThreshold':self.mutthreshold,
                    'QAverage':self.obj_cfgparser.get('fastq_filter_translate_ref', 'qaverage'),
                    'QLimit':self.obj_cfgparser.get('fastq_filter_translate_ref', 'qlimit'),
                    'In_File':file_input_fastq_ref,
                    'Out_Prefix':self.directory + self.output_prefix + "_Ref",
                    'Enable_Anchors':self.obj_cfgparser.get('fastq_filter_translate_ref', 'enable_anchors')
                    }
                except NoSectionError:
                    print("[Protocols:Fitness Error] The fastq_filter_translate_ref config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [fastq_filter_translate_ref] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The fastq_filter_translate_ref config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()

                #Create the object then call the first file
                obj_fastqfilter_ref = fastq_filter_translate(dict_fastqread_options_ref)

                print("[Protocols:Fitness] Reading the reference fastq files.")
                file_output.write("[Protocols:Fitness] Reading the reference fastq files.\n")
                file_output.write(obj_fastqfilter_ref.fastq_filter_translate() + "\n")

            if self.dict_workflow['fastq_filter_translate_sel']:

                try:
                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('fastq_filter_translate_sel', 'fastq_file')) > 0:
                        file_input_fastq_sel = self.obj_cfgparser.get('fastq_filter_translate_sel', 'fastq_file')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Sel_Merge.fastq"):
                            file_input_fastq_sel = self.directory + self.output_prefix + "_Sel_Merge.fastq"

                    dict_fastqread_options_sel = {
                    'WTDNA':self.wtdna,
                    'WTAA':self.wtaa,
                    'FirstAAMutated':self.firstaamutated,
                    'LastAAMutated':self.lastaamutated,
                    'Processes':self.processes,
                    '5pAnchor':self.obj_cfgparser.get('fastq_filter_translate_sel', 'fiveprimeanchor'),
                    'MutThreshold':self.mutthreshold,
                    'QAverage':self.obj_cfgparser.get('fastq_filter_translate_sel', 'qaverage'),
                    'QLimit':self.obj_cfgparser.get('fastq_filter_translate_sel', 'qlimit'),
                    'In_File':file_input_fastq_sel,
                    'Out_Prefix':self.directory + self.output_prefix + "_Sel",
                    'Enable_Anchors':self.obj_cfgparser.get('fastq_filter_translate_sel', 'enable_anchors')
                    }

                except NoSectionError:
                    print("[Protocols:Fitness Error] The fastq_filter_translate_sel config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [fastq_filter_translate_sel] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The fastq_filter_translate_sel config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()      

                #Create the object then call the first file
                obj_fastqfilter_sel = fastq_filter_translate(dict_fastqread_options_sel)

                print("[Protocols:Fitness] Reading the selected fastq files.")
                file_output.write("[Protocols:Fitness] Reading the selected fastq files.\n")
                file_output.write(obj_fastqfilter_sel.fastq_filter_translate() + "\n")

            """
            *****************************************
            Filter Counter Section
            *****************************************
            """            
            if self.dict_workflow['filter_counter_sel'] or self.dict_workflow['filter_counter_ref']:
                #Import our filter_counter file
                try:
                    from pact.sequencing.mut_filter import mut_filter
                    from pact.sequencing.mut_counter import mut_counter
                except ImportError:
                    print("[Protocols:Fitness Error] mut_filter was not found.")        

            if self.dict_workflow['filter_counter_ref']:
                try:
                    if len(self.obj_cfgparser.get('filter_counter_ref', 'read_file')) > 0:
                        file_input_fastqread_ref = self.obj_cfgparser.get('filter_counter_ref', 'read_file')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Ref_Read.tsv"):
                            file_input_fastqread_ref = self.directory + self.output_prefix + "_Ref_Read.tsv"

                    dict_filtercounter_options_ref = {
                    'WTDNA':self.wtdna,
                    'WTAA':self.wtaa,
                    'FirstAAMutated':self.firstaamutated,
                    'LastAAMutated':self.lastaamutated,
                    'MutThreshold':self.mutthreshold,
                    'Processes':self.processes,
                    'In_File':file_input_fastqread_ref,
                    'Out_Prefix':self.directory + self.output_prefix + "_Ref",
                    'mutcodons':self.mutcodons,
                    }
                except NoSectionError:
                    print("[Protocols:Fitness Error] The filter_counter config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [filter_counter_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The filter_counter config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()

                #Create the object then call the first file
                obj_filter_ref = mut_filter(dict_filtercounter_options_ref)

                print("[Protocols:Fitness] Filtering the reference tsv.")
                file_output.write("[Protocols:Fitness] Filtering the reference tsv.\n")
                file_output.write(obj_filter_ref.mut_filter() + "\n")

                #Create the object then call the first file
                obj_counter_ref = mut_counter(dict_filtercounter_options_ref)

                print("[Protocols:Fitness] Counting the reference mutations.")
                file_output.write("[Protocols:Fitness] Counting the reference mutations.\n")
                file_output.write(obj_counter_ref.mut_counter() + "\n")

            if self.dict_workflow['filter_counter_sel']:

                try:
                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('filter_counter_sel', 'read_file')) > 0:
                        file_input_fastqread_sel = self.obj_cfgparser.get('filter_counter_sel', 'read_file')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Sel_Read.tsv"):
                            file_input_fastqread_sel = self.directory + self.output_prefix + "_Sel_Read.tsv"

                    dict_filtercounter_options_sel = {
                    'WTDNA':self.wtdna,
                    'WTAA':self.wtaa,
                    'FirstAAMutated':self.firstaamutated,
                    'LastAAMutated':self.lastaamutated,
                    'MutThreshold':self.mutthreshold,
                    'Processes':self.processes,
                    'In_File':file_input_fastqread_sel,
                    'Out_Prefix':self.directory + self.output_prefix + "_Sel",
                    'mutcodons':self.mutcodons,
                    }

                except NoSectionError:
                    print("[Protocols:Fitness Error] The filter_counter config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [filter_counter_X] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The filter_counter config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()

                #Create the object then call the first file
                obj_filter = mut_filter(dict_filtercounter_options_sel)

                print("[Protocols:Fitness] Filtering the selected tsv.")
                file_output.write("[Protocols:Fitness] Filtering the selected tsv.\n")
                file_output.write(obj_filter.mut_filter() + "\n")

                #Create the object then call the first file
                obj_counter_sel = mut_counter(dict_filtercounter_options_sel)

                print("[Protocols:Fitness] Counting the selected mutations.")
                file_output.write("[Protocols:Fitness] Counting the selected mutations.\n")
                file_output.write(obj_counter_sel.mut_counter() + "\n")
        
            """
            *****************************************
            Enrichment Section
            *****************************************
            """
            if self.dict_workflow['enrichment']:

                # do we consider rejected mutations (based on design) from the total count?
                try:
                    consider_rejected = self.obj_cfgparser.get('enrichment', 'consider_rejected')
                except NoOptionError:
                    consider_rejected = 'false'
                    
                try:
                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('enrichment', 'ref_count_wildtype')) > 0:
                        file_input_countwt_ref = self.obj_cfgparser.get('enrichment', 'ref_count_wildtype')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Ref_Counted_WildType.tsv"):
                            file_input_countwt_ref = self.directory + self.output_prefix + "_Ref_Counted_WildType.tsv"

                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('enrichment', 'sel_count_wildtype')) > 0:
                        file_input_countwt_sel = self.obj_cfgparser.get('enrichment', 'sel_count_wildtype')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Sel_Counted_WildType.tsv"):
                            file_input_countwt_sel = self.directory + self.output_prefix + "_Sel_Counted_WildType.tsv"


                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('enrichment', 'ref_count')) > 0:
                        file_input_count_ref = self.obj_cfgparser.get('enrichment', 'ref_count')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Ref_Counted.tsv"):
                            file_input_count_ref = self.directory + self.output_prefix + "_Ref_Counted.tsv"

                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('enrichment', 'sel_count')) > 0:
                        file_input_count_sel = self.obj_cfgparser.get('enrichment', 'sel_count')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Sel_Counted.tsv"):
                            file_input_count_sel = self.directory + self.output_prefix + "_Sel_Counted.tsv"


                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('enrichment', 'ref_count_rejected')) > 0:
                        file_input_countrej_ref = self.obj_cfgparser.get('enrichment', 'ref_count_rejected')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Ref_Counted_Rejected.tsv"):
                            file_input_countrej_ref = self.directory + self.output_prefix + "_Ref_Counted_Rejected.tsv"

                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('enrichment', 'sel_count_rejected')) > 0:
                        file_input_countrej_sel = self.obj_cfgparser.get('enrichment', 'sel_count_rejected')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_Sel_Counted_Rejected.tsv"):
                            file_input_countrej_sel = self.directory + self.output_prefix + "_Sel_Counted_Rejected.tsv"

                    enrichment_settings = {
                    'WTDNA':self.wtdna,
                    'FirstAAMutated':self.firstaamutated,
                    'LastAAMutated':self.lastaamutated,
                    'Out_Prefix':self.directory + self.output_prefix,
                    'mutcodons':self.mutcodons,

                    "Ref_Count_WildType":file_input_countwt_ref,
                    "Sel_Count_WildType":file_input_countwt_sel,
                    "Ref_Count":file_input_count_ref,
                    "Sel_Count":file_input_count_sel,
                    "Ref_Count_Rejected":file_input_countrej_ref,
                    "Sel_Count_Rejected":file_input_countrej_sel,

                    "Ref_Count_Threshold":self.obj_cfgparser.get('enrichment', 'ref_count_threshold'),
                    "Sel_Count_Threshold":self.obj_cfgparser.get('enrichment', 'sel_count_threshold'),
                    "Strict_Count_Threshold":self.obj_cfgparser.get('enrichment', 'strict_count_threshold'),
                    "consider_rejected":consider_rejected,
                    }
                except NoSectionError:
                    print("[Protocols:Fitness Error] The enrichment config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [enrichment] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The enrichment config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()

                #Import our enrichment_fitness file
                try:
                    from pact.sequencing.enrichment import enrichment
                except ImportError:
                    print("[Protocols:Fitness Error] enrichment was not found.")

                #Create the object
                obj_enrichment = enrichment(enrichment_settings)
    
                #Calculate the enrichment
                print("[Protocols:Fitness] Calculating the enrichment")
                file_output.write("[Protocols:Fitness] Calculating the enrichment\n" + obj_enrichment.enrichment())

            """
            *****************************************
            Fitness Section
            *****************************************
            """
            if self.dict_workflow['fitness']:

                #Import our fitness file
                try:
                    from pact.sequencing.fitness import fitness
                except ImportError:
                    print("[Protocols:Fitness Error] fitness was not found.")

                #Create the object
                obj_fitness = fitness(self.obj_cfgparser, self.dict_programs, 
                                       {'directory':self.directory,
                                        'WTDNA':self.wtdna,
                                        'WTAA':self.wtaa,
                                        'FirstAAMutated':self.firstaamutated,
                                        'LastAAMutated':self.lastaamutated,
                                        'Out_Prefix':self.directory + self.output_prefix,
                                        'mutcodons':self.mutcodons,
                                        'library_type':self.mutationtype})
    
                #Calculate the fitness
                print("[Protocols:Fitness] Calculating the fitness")
                file_output.write("[Protocols:Fitness] Calculating the fitness\n" + obj_fitness.fitness())

            """
            *****************************************
            Calculate mutation freqs and mutual info
            *****************************************
            """
            if self.mutationtype == 'multiple':
                if self.dict_workflow['multiple_freq_mi']:
                    #Import our codon frequency and mutual information class
                    try:
                        from pact.sequencing.multiple_freq_mi import multiple_freq_mi
                    except ImportError:
                        print("[Protocols:Fitness Error] multiple_freq_mi was not found.")

                    #Create the object
                    obj_freqmi = multiple_freq_mi(self.obj_cfgparser, self.dict_programs, 
                                           {'directory':self.directory,
                                            'WTDNA':self.wtdna,
                                            'WTAA':self.wtaa,
                                            'FirstAAMutated':self.firstaamutated,
                                            'LastAAMutated':self.lastaamutated,
                                            'Out_Prefix':self.directory + self.output_prefix,
                                            'mutcodons':self.mutcodons,
                                            'library_type':self.mutationtype})
    
                    #Calculate the codon frequency and mutual information
                    print("[Protocols:Fitness] Calculating the codon frequency and mutual information")
                    file_output.write("[Protocols:Fitness] Calculating the fitness\n" + obj_freqmi.multiple_freq_mi())

            """
            *****************************************
            Library Stats Section
            *****************************************
            """
            if self.dict_workflow['library_stats']:
            
                try:
                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('library_stats', 'pact_enrichment_summary')) > 0:
                        file_enrich_summary = self.obj_cfgparser.get('library_stats', 'pact_enrichment_summary')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_enrichment_summary.pact"):
                            file_enrich_summary = self.directory + self.output_prefix + "_enrichment_summary.pact"

                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('library_stats', 'pact_fitness_nonsynon')) > 0:
                        file_fitness_nonsynon = self.obj_cfgparser.get('library_stats', 'pact_fitness_nonsynon')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_fitness_nonsynon.pact"):
                            file_fitness_nonsynon = self.directory + self.output_prefix + "_fitness_nonsynon.pact"

                    #Test if we have a special input file
                    if len(self.obj_cfgparser.get('library_stats', 'pact_fitness_wtsynon')) > 0:
                        file_fitness_wtsynon = self.obj_cfgparser.get('library_stats', 'pact_fitness_wtsynon')
                    else:
                        if file_checker(self.directory + self.output_prefix + "_fitness_wtsynon.pact"):
                            file_fitness_wtsynon = self.directory + self.output_prefix + "_fitness_wtsynon.pact"

                    library_settings = {
                        'file_summary':file_enrich_summary,
                        'pact_fitness_nonsynon':file_fitness_nonsynon,
                        'pact_fitness_wtsynon':file_fitness_wtsynon,

                        'WTDNA':self.wtdna,
                        'WTAA':self.wtaa,
                        'FirstAAMutated':self.firstaamutated,
                        'LastAAMutated':self.lastaamutated,
                        'Out_Prefix':self.directory + self.output_prefix,
                        'mutcodons':self.mutcodons,
                        'library_type':self.mutationtype,

                        'codon_type':self.obj_cfgparser.get('library_stats', 'codon_type'),
                        }
                except NoSectionError:
                    print("[Protocols:Fitness Error] The library_stats config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the [library_stats] section.")
                    quit()
                except NoOptionError:
                    print("[Protocols:Fitness Error] The library_stats config file is incorrect.")
                    print("[Protocols:Fitness Error] There is something wrong with the name of a option flag.")
                    quit()

                #Import our library_stats file
                try:
                    from pact.sequencing.library_stats import library_stats
                except ImportError:
                    print("[Protocols:Fitness Error] library_stats was not found.")

                #Create or library stats object
                obj_libstat = library_stats(library_settings)
    
                #Call our entrypoint method
                file_output.write(obj_libstat.library_stats() + "\n")

        return

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Protocols:" + str_protocol_name + " Error] This protocol needs to be ran within the context of PACT.")
