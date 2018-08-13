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

"""residue_chemical_size - Return the amino acid classifier of chemical and size for a mutation pair"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Residue Chemical Size Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from pact.pact_common import pretty_counter_dicts

#Set the author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

class residue_chemical_size:
    """Return the amino acid classifier of chemical and size for a mutation pair"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Setup the globals"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        return

    def chemical_polarity(self, from_aa, to_aa):
        """
        Polar/Charged/Non-Polar
        CR = Charge Reversal
        CC = Charged to Charged
        CP = Charged to Polar
        CN = Charged to Non-Polar
        PC = Polar to Charged
        PP = Polar to Polar
        PN = Polar to Non-Polar
        NC = Non-Polar to Charged
        NP = Non-Polar to Polar
        NN = Non-Polar to Non-Polar
        Charge Neg   D/E
        Charge Pos   R/K/H
        Polar        S/T/N/Q/C
        Non-Polar    G/A/L/I/M/V/F/W/Y/P
        * = * = Stop
        """

        #Setup the transform dict
        dict_type = {
            "D":"-", "E":"-", "R":"+", "K":"+", "H":"+",
            "S":"P", "T":"P", "N":"P", "Q":"P", "C":"P",
            "G":"N", "A":"N", "L":"N", "I":"N", "M":"N",
            "V":"N", "F":"N", "W":"N", "Y":"N", "P":"N",
            "*":"*"}
        
        #Setup the return lookup dict
        dict_polarity = {
            "-+":"CR", "+-":"CR", "--":"CC", "++":"CC",
            "-P":"CP", "+P":"CP", "-N":"CN", "+N":"CN",
            "P-":"PC", "P+":"PC", "PP":"PP", "PN":"PN",
            "N-":"NC", "N+":"NC", "NN":"NN", "NP":"NP",
            "*":"*"}

        #Return our classifier code
        return dict_polarity[dict_type[from_aa.upper()] + dict_type[to_aa.upper()]]

    def chemical_aromatics(self, from_aa, to_aa):
        """
        Aromatic/Non-Aromatic
        AN = Aromatic to Non-Aromatic
        NA = Non-Aromatic to Aromatic
        AA = Aromatic to Aromatic
        NN = Non-Aromatic to Non_Aromatic
        Aromatic     F/W/Y/H
        Non_aromatic D/E/R/K/S/T/N/Q/C/G/A/L/I/M/V/P
        """

        #Setup the transform dict
        dict_type = {
            "D":"N", "E":"N", "R":"N", "K":"N", "H":"A",
            "S":"N", "T":"N", "N":"N", "Q":"N", "C":"N",
            "G":"N", "A":"N", "L":"N", "I":"N", "M":"N",
            "V":"N", "F":"A", "W":"A", "Y":"A", "P":"N"}

        #Return our classifier code
        return dict_type[from_aa.upper()] + dict_type[to_aa.upper()]

    def chemical_hydro_phob_phil(self, from_aa, to_aa):
        """
        Hydrophobic/Hydrophilic/Cysteine/Proline
        CC = Cys to Cys
        PP = Pro to Pro
        CP = Cys to Pro
        PC = Pro to Cys
        BL = hydrophoBic to hydrophiLic
        LB = hydrophiLic to hydrophoBic
        BB = hydrophoBic
        LL = hydrophiLic
        
        CL = Cys to hydrophiLic
        CB = Cys to hydrophoBic
        LC = hydrophiLic to Cys
        BC = hydrophoBic to Cys

        PL = Pro to hydrophiLic
        PB = Pro to hydrophoBic
        LP = hydrophiLic to Pro
        BP = hydrophoBic to Pro
        Hydrophobic  A/L/I/M/V/F/W/Y/G
        Hydrophilic  D/E/R/K/H/S/T/N/Q
        Cysteine     C
        Proline      P
        """

        #Setup the transform dict
        dict_type = {
            "D":"L", "E":"L", "R":"L", "K":"L", "H":"L",
            "S":"L", "T":"L", "N":"L", "Q":"L", "C":"C",
            "G":"B", "A":"B", "L":"B", "I":"B", "M":"B",
            "V":"B", "F":"B", "W":"B", "Y":"B", "P":"P",
            "*":"*"}

        #Return our classifier code
        return dict_type[from_aa.upper()] + dict_type[to_aa.upper()]

    def size_change(self, from_aa, to_aa):
        """
        Size
        SB = Small to Big
        BS = Big to Small
        SS = Small to Small
        BB = Big to Big
        PX = PROLINE to X
        XP = X to PROLINE

        Small   D/S/T/V/N/C/A/G
        Big     R/K/E/L/I/M/Q/F/W/Y/H
        Pro     P
        """

        #Setup the transform dict
        dict_type = {
            "D":"S", "E":"B", "R":"B", "K":"B", "H":"B",
            "S":"S", "T":"S", "N":"S", "Q":"B", "C":"S",
            "G":"S", "A":"S", "L":"B", "I":"B", "M":"B",
            "V":"S", "F":"B", "W":"B", "Y":"B", "P":"P"}

        #Return our classifier code
        return dict_type[from_aa.upper()] + dict_type[to_aa.upper()]

    def hydropathy_change(self, from_aa, to_aa):
        """
        Calculate the change in the hydropathy index (Kyte 1982 JMB)
        Postive = Hydrophobic, Negative = Hydrophilic
        Modify with Wimley White 1006 Nature SMB
        """

        #Setup the transform dict
        dict_type = {
            "D":-3.5, "E":-3.5, "R":-4.5, "K":-3.9, "H":-3.2,
            "S":-0.8, "T":-0.7, "N":-3.5, "Q":-3.5, "C":2.5,
            "G":-0.4, "A":1.8, "L":3.8, "I":4.5, "M":1.9,
            "V":4.2, "F":2.8, "W":-0.9, "Y":-1.3, "P":-1.6}

        #Return our classifier code
        return round(dict_type[to_aa.upper()] - dict_type[from_aa.upper()], 2)

    def process_dataset(self, dict_merged_datasets):
        """Process the fitness dataset and return the location, mutation, data"""

        dict_return = {}

        #Loop all of the datasets
        for dataset in dict_merged_datasets:

            dict_return[dataset] = {}

            #Loop the locations
            for key in dict_merged_datasets[dataset]:

                dict_return[dataset][key] = {}

                #Loop the muts
                for mut in dict_merged_datasets[dataset][key]:

                    dict_return[dataset][key][mut] = {}

                    wt_residue = dict_merged_datasets[dataset][key][mut]['wt_residue']

                    #Handle WT
                    if mut == wt_residue:
                        dict_return[dataset][key][mut] = {'polarity':'WT', 'aromatics':'WT',
                                                 'philic_phobic':'WT', 'size':'WT',
                                                 'hydropathy':'WT'}
                        continue

                    #Handle Stop
                    if mut == "*":
                        dict_return[dataset][key][mut] = {'polarity':'STOP', 'aromatics':'STOP',
                                                 'philic_phobic':'STOP', 'size':'STOP',
                                                 'hydropathy':'STOP'}
                        continue

                    dict_return[dataset][key][mut] = {
                        'polarity':self.chemical_polarity(wt_residue, mut),
                        'aromatics':self.chemical_aromatics(wt_residue, mut),
                        'philic_phobic':self.chemical_hydro_phob_phil(wt_residue, mut),
                        'size':self.size_change(wt_residue, mut),
                        'hydropathy':self.hydropathy_change(wt_residue, mut)
                        }

        return dict_return

    def mut_info(self, from_aa, to_aa):
        """Return the entire info on a per mutation basis"""
        return {
                'polarity':self.chemical_polarity(from_aa, to_aa),
                'aromatics':self.chemical_aromatics(from_aa, to_aa),
                'philic_phobic':self.chemical_hydro_phob_phil(from_aa, to_aa),
                'size':self.size_change(from_aa, to_aa),
                'hydropathy':self.hydropathy_change(from_aa, to_aa)
                }

    def classified_count(self, dict_rcs, dataset, dict_classified):
        """Count our classifiers"""

        #Import our counter
        from collections import Counter

        """
        Polar/Charged to Polar/Charged
        Charge Reversal
        Polar/Charge to Hydrophobic/Aromatic
        Hydrophobic/Aromatic to Polar/Charged
        To/From Proline
        Hydrophobic/Aromatic to Hydrophobic/Aromatic

        Big to Big
        Big to Small
        To/From Proline
        Small to Big
        Small to Small
        """

        #Setup our lists
        list_pctopc = []
        list_cr = []
        list_pctoha = []
        list_hatopc = []
        list_tofrompro = []
        list_hatoha = []
        list_bb = []
        list_bs = []
        list_tofrompro_size = []
        list_sb = []
        list_ss = []

        #Loop the locations
        for loc in dict_classified:

            #Loop the mutations
            for mut in dict_classified[loc]:

                #Polar/Charged to Polar/Charged 'philic_phobic'[LL]
                if dict_rcs[dataset][loc][mut]['philic_phobic'] == "LL":
                    list_pctopc.append(dict_classified[loc][mut])

                #Charge Reversal 'polarity'[CR]
                if dict_rcs[dataset][loc][mut]['polarity'] == "CR":
                    list_cr.append(dict_classified[loc][mut])

                #Polar/Charge to Hydrophobic/Aromatic 'philic_phobic'[LB, LC]
                if (dict_rcs[dataset][loc][mut]['philic_phobic'] == "LB" or
                    dict_rcs[dataset][loc][mut]['philic_phobic'] == "LC"):
                    list_pctoha.append(dict_classified[loc][mut])

                #Hydrophobic/Aromatic to Polar/Charged 'philic_phobic'[BL, CL]
                if (dict_rcs[dataset][loc][mut]['philic_phobic'] == "BL" or
                    dict_rcs[dataset][loc][mut]['philic_phobic'] == "CL"):
                    list_hatopc.append(dict_classified[loc][mut])

                #To/From Proline 'philic_phobic'[PP, CP, PC, PL, PB, LP, BP]
                if "P" in dict_rcs[dataset][loc][mut]['philic_phobic']:
                    list_tofrompro.append(dict_classified[loc][mut])

                #Hydrophobic/Aromatic to Hydrophobic/Aromatic 'philic_phobic'[BB, CC, BC, CB]
                if (dict_rcs[dataset][loc][mut]['philic_phobic'] == "BB" or
                    dict_rcs[dataset][loc][mut]['philic_phobic'] == "CC" or
                    dict_rcs[dataset][loc][mut]['philic_phobic'] == "BC" or
                    dict_rcs[dataset][loc][mut]['philic_phobic'] == "CB"):
                    list_hatoha.append(dict_classified[loc][mut])

                #Big to Big 'size'[BB]
                if dict_rcs[dataset][loc][mut]['size'] == "BB":
                    list_bb.append(dict_classified[loc][mut])

                #Big to Small 'size'[BS]
                if dict_rcs[dataset][loc][mut]['size'] == "BS":
                    list_bs.append(dict_classified[loc][mut])

                #To/From Proline 'size'[PB, PS, SP, BP]
                if "P" in dict_rcs[dataset][loc][mut]['size']:
                    list_tofrompro_size.append(dict_classified[loc][mut])

                #Small to Big 'size'[SB]
                if dict_rcs[dataset][loc][mut]['size'] == "SB":
                    list_sb.append(dict_classified[loc][mut])

                #Small to Small 'size'[SS]
                if dict_rcs[dataset][loc][mut]['size'] == "SS":
                    list_ss.append(dict_classified[loc][mut])

        #Count the lists
        str_return = '\n'.join(map(str, [
            "****Chemical Change****",
            "Polar/Charged to Polar/Charged",
            pretty_counter_dicts(dict(Counter(list_pctopc))),
            "",
            "Charge Reversal",
            pretty_counter_dicts(dict(Counter(list_cr))),
            "",
            "Polar/Charge to Hydrophobic/Aromatic",
            pretty_counter_dicts(dict(Counter(list_pctoha))),
            "",
            "Hydrophobic/Aromatic to Polar/Charged",
            pretty_counter_dicts(dict(Counter(list_hatopc))),
            "",
            "To/From Proline",
            pretty_counter_dicts(dict(Counter(list_tofrompro))),
            "",
            "Hydrophobic/Aromatic to Hydrophobic/Aromatic",
            pretty_counter_dicts(dict(Counter(list_hatoha))),
            "",
            "****Size Change****",
            "Big to Big",
            pretty_counter_dicts(dict(Counter(list_bb))),
            "",
            "Big to Small",
            pretty_counter_dicts(dict(Counter(list_bs))),
            "",
            "To/From Proline",
            pretty_counter_dicts(dict(Counter(list_tofrompro_size))),
            "",
            "Small to Big",
            pretty_counter_dicts(dict(Counter(list_sb))),
            "",
            "Small to Small",
            pretty_counter_dicts(dict(Counter(list_ss))),
            ]))

        print(str_return)
        return str_return

if __name__ == '__main__':
    #Remind the user that the classifier needs to be ran within the context of PACT
    print("[Residue Chemical Size Error] This classifier needs to be ran within the context of PACT.")   
