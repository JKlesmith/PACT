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

"""dataset vs dataset - analyze and plot two pact datasets against each other"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[Set vs Set Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from configparser import NoOptionError
from pact.pact_common import file_checker
from scipy.stats import linregress
from statistics import mean, stdev

#The author information
__author__ = "Justin R. Klesmith"
__copyright__ = "Copyright 2018, Justin R. Klesmith"
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "justinklesmith@gmail.com"]

class set_vs_set:
    """dataset vs dataset - analyze and plot two pact datasets against each other"""

    def __init__(self, obj_config, dict_programs, dict_protocolconfig):
        """Initialize the class varibles"""
        
        #Get the config file parser object from the calling class
        self.config_file = obj_config

        #Get the PACT dict programs
        self.dict_programs = dict_programs

        #Get the extra config options assigned by the protocol
        self.dict_protocolconfig = dict_protocolconfig

        #Set the output prefix
        self.directory = self.dict_protocolconfig['directory']
        self.output_prefix = self.config_file.get('global', 'output_prefix')

        #Classifier specific variables
        self.aa_table = '*FWYPMILVAGCSTNQDEHKR'

        #Options
        try:
            self.dataset_x = self.config_file.get('setvsset', 'dataset_x')
            self.dataset_y = self.config_file.get('setvsset', 'dataset_y')
            
            self.x_column = self.config_file.get('setvsset', 'x_column').lower()
            self.y_column = self.config_file.get('setvsset', 'y_column').lower()

            self.x_label = self.config_file.get('setvsset', 'x_axis_label')
            self.y_label = self.config_file.get('setvsset', 'y_axis_label')

            self.ref_threshold = int(self.config_file.get('setvsset', 'ref_threshold'))
            self.sel_threshold = int(self.config_file.get('setvsset', 'sel_threshold'))

            self.x_min = float(self.config_file.get('setvsset', 'x_axis_min'))
            self.x_max = float(self.config_file.get('setvsset', 'x_axis_max'))
            self.y_min = float(self.config_file.get('setvsset', 'y_axis_min'))
            self.y_max = float(self.config_file.get('setvsset', 'y_axis_max'))
            
            self.xy_type = self.config_file.get('setvsset', 'xy_scatter').lower()
            self.xy_grouping = self.config_file.get('setvsset', 'xy_scatter_type').lower()
            self.output_csv = self.config_file.get('setvsset', 'output_csv').lower()

            self.pact_headless = self.config_file.get('setvsset', 'headless').lower()

            self.outlier_threshold = float(self.config_file.get('setvsset', 'outlier_threshold'))
            self.winner_threshold = float(self.config_file.get('setvsset', 'winner_threshold'))
            self.amino_list = self.config_file.get('setvsset', 'amino_acid_highlight').upper().split(',')
            self.point_color = self.config_file.get('setvsset', 'point_color').lower()

            if float(self.config_file.get('setvsset', 'sd_boundaries')) > 0:
                self.sd_bounds = float(self.config_file.get('setvsset', 'sd_boundaries'))
            else:
                self.sd_bounds = None

            #Set the regression option
            self.bool_regression = self.config_file.get('setvsset', 'regression').lower()
            self.onetoone = self.config_file.get('setvsset', '1to1line').lower()
            self.count_shared = self.config_file.get('setvsset', 'shared_counts').lower()

        except NoOptionError:
            print("[Set vs Set Error] Missing config file elements.")
            quit()

        #Check the validity of the groups
        for letter in self.amino_list:
            if letter not in self.aa_table:
                print("[Set vs Set Error] Unknown letter in amino_acid_highlight.")

        return

    def xy_list(self, dict_merged_datasets):
        """Return two lists of x and y"""

        #Make a blank list
        x = []
        y = []
        x_nan = []
        y_nan = []

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    x_nan.append(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]))
                    y_nan.append(float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Add our mutation
                x.append(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]))
                y.append(float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))

        return x, y, x_nan, y_nan

    def xy_list_grouplocout(self, dict_merged_datasets, dict_color = {}):
        """Return two lists of x and y"""

        #Make a blank list
        x = []
        y = []
        x_nan = []
        y_nan = []
        xerr = []
        yerr = []
        color = []
        size = []

        #Print a output table
        str_outlier = ','.join(map(str, [
            "Difference:",
            "Calculated_Diff",
            "Location",
            "Wild-Type",
            "Mutation",
            "X Fitness",
            "X Column",
            "Y Fitness",
            "Y Column"
            ])) + "\n"

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Setup the temp lists to add to
            tempx = []
            tempy = []

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip Stop
                if mut == "*":
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    x_nan.append(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]))
                    y_nan.append(float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Add our mutation
                unity_difference = abs(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]) -
                float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))

                #Report if our difference is greater than our threshold (1=1 line)
                if unity_difference > self.outlier_threshold:
                    tempx.append(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]))
                    tempy.append(float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))

            #Calculate the average and SD
            if len(tempx) == 0 or len(tempy) == 0:
                continue

            x.append(mean(tempx))
            y.append(mean(tempy))

            if len(tempx) < 2:
                xerr.append(0)
            else:
                xerr.append(stdev(tempx))

            if len(tempy) < 2:
                yerr.append(0)
            else:
                yerr.append(stdev(tempy))

            #Parse our dict
            if (key not in dict_color or mut not in dict_color[key]):
                color.append('black')
            else:
                color.append(dict_color[key][mut])

            #Calculate our size
            size.append(2**(len(tempx) * 0.5))

            #Make a report
            str_outlier = str_outlier + ','.join(map(str, [
                "Difference:",
                key,
                round(x[-1], 3),
                round(xerr[-1], 3),
                round(y[-1], 3),
                round(yerr[-1], 3),
                ])) + "\n"

        print(str_outlier)

        return x, y, x_nan, y_nan, xerr, yerr, color, size


    def dictcolor_to_list(self, dict_merged_datasets, dict_color):
        """Convert a custom color dict to a list"""

        #Make a list for colors
        list_color = []

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Parse our dict
                if (key not in dict_color or mut not in dict_color[key]):
                    list_color.append('black')
                else:
                    list_color.append(dict_color[key][mut])

        return list_color

    def dictsize_to_list(self, dict_merged_datasets, dict_size):
        """Convert a custom size dict to a list"""

        #Make a list for colors
        list_size = []

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Parse our dict
                if (key not in dict_size or mut not in dict_size[key]):
                    list_size.append(0)
                else:
                    list_size.append(dict_size[key][mut])

        return list_size

    def dicterror_to_list(self, dict_merged_datasets, dict_error):
        """Convert a custom size dict to a list"""

        #Make a list for colors
        list_xerr = []
        list_yerr = []

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Parse our dict
                if (key not in dict_error or mut not in dict_error[key]):
                    list_xerr.append(0)
                    list_yerr.append(0)
                else:
                    list_xerr.append(dict_error[key][mut]['xerr'])
                    list_yerr.append(dict_error[key][mut]['yerr'])

        return list_xerr, list_yerr


    def outliers(self, dict_merged_datasets):
        """Return a list of points that differ by a significant amount vs a fit"""

        #Make a list for colors
        list_color = []

        #Print a output table
        str_outlier = ','.join(map(str, [
            "Difference:",
            "Calculated_Diff",
            "Location",
            "Wild-Type",
            "Mutation",
            "X Fitness",
            "X Column",
            "Y Fitness",
            "Y Column"
            ])) + "\n"

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Perform our check
                unity_difference = abs(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]) -
                float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))

                #Report if our difference is greater than our threshold (1=1 line)
                if unity_difference > self.outlier_threshold:

                    str_outlier = str_outlier + ','.join(map(str, [
                        "Difference:",
                        round(unity_difference, 3),
                        key,
                        dict_merged_datasets[self.dataset_x][key][mut]['wt_residue'],
                        mut,
                        round(dict_merged_datasets[self.dataset_x][key][mut]['fitness'], 3),
                        round(dict_merged_datasets[self.dataset_x][key][mut][self.x_column], 3),
                        round(dict_merged_datasets[self.dataset_y][key][mut]['fitness'], 3),
                        round(dict_merged_datasets[self.dataset_y][key][mut][self.y_column], 3),
                        ])) + "\n"

                    list_color.append("red")
                else:
                    list_color.append("black")

        print(str_outlier)

        return list_color, str_outlier

    def outlier_sign(self, dict_merged_datasets):
        """Return a list of points that differ by a significant amount vs a fit"""

        #Make a list for colors
        list_color = []

        #Print a output table
        str_outlier = ','.join(map(str, [
            "Difference:",
            "Calculated_Diff",
            "Location",
            "Wild-Type",
            "Mutation",
            "X Fitness",
            "X Column",
            "Y Fitness",
            "Y Column"
            ])) + "\n"

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Perform our check
                unity_difference = abs(float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]) -
                float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]))

                #Report if our difference is greater than our threshold (1=1 line)
                if unity_difference > self.outlier_threshold:

                    str_outlier = str_outlier + ','.join(map(str, [
                        "Difference:",
                        round(unity_difference, 3),
                        key,
                        dict_merged_datasets[self.dataset_x][key][mut]['wt_residue'],
                        mut,
                        round(dict_merged_datasets[self.dataset_x][key][mut]['fitness'], 3),
                        round(dict_merged_datasets[self.dataset_x][key][mut][self.x_column], 3),
                        round(dict_merged_datasets[self.dataset_y][key][mut]['fitness'], 3),
                        round(dict_merged_datasets[self.dataset_y][key][mut][self.y_column], 3),
                        ])) + "\n"

                    #Check for the sign difference
                    if ((dict_merged_datasets[self.dataset_x][key][mut][self.x_column] > 0 and 
                        dict_merged_datasets[self.dataset_y][key][mut][self.y_column] > 0) or 
                        (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] < 0 and 
                        dict_merged_datasets[self.dataset_y][key][mut][self.y_column] < 0)):
                        list_color.append("red")
                    else:
                        list_color.append("blue")
                else:
                    list_color.append("black")

        print(str_outlier)

        return list_color, str_outlier

    def winners(self, dict_merged_datasets):
        """Return a list of winning points above a certain threshold"""

        #Make a list for colors
        list_color = []

        #Print a output table
        str_winner = ','.join(map(str, [
            "Winner:",
            "Location",
            "Wild-Type",
            "Mutation",
            "X Fitness",
            "X Column",
            "Y Fitness",
            "Y Column"
            ])) + "\n"

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Report if our difference is greater than our threshold (1=1 line)
                if (float(dict_merged_datasets[self.dataset_x][key][mut][self.x_column]) >= self.winner_threshold and 
                     float(dict_merged_datasets[self.dataset_y][key][mut][self.y_column]) >= self.winner_threshold):

                    str_winner = str_winner + ','.join(map(str, [
                        "Winner:",
                        key,
                        dict_merged_datasets[self.dataset_x][key][mut]['wt_residue'],
                        mut,
                        round(dict_merged_datasets[self.dataset_x][key][mut]['fitness'], 3),
                        round(dict_merged_datasets[self.dataset_x][key][mut][self.x_column], 3),
                        round(dict_merged_datasets[self.dataset_y][key][mut]['fitness'], 3),
                        round(dict_merged_datasets[self.dataset_y][key][mut][self.y_column], 3),
                        ])) + "\n"

                    list_color.append("red")
                else:
                    list_color.append("black")

        print(str_winner)

        return list_color, str_winner

    def aminos(self, dict_merged_datasets):
        """Return a list of amino acid points"""

        #Make a list for colors
        list_color = []

        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Report if our difference is greater than our threshold (1=1 line)
                if mut in self.amino_list:
                    list_color.append("red")
                else:
                    list_color.append("black")

        return list_color


    def xy_scatter_error(self, x, y, color = [], size = [], xerr = [], yerr = []):
        """Plot the scatter with errorbars"""

        #Import the graphing toolkit
        from matplotlib import use, rcParams

        #Set the matplotlib base settings (headless, font)
        if self.pact_headless == "true":
            use('Agg')

        rcParams['font.sans-serif']='Arial'
        rcParams['font.size']='22'

        from matplotlib import pyplot as plt
        from matplotlib.ticker import FormatStrFormatter

        #Check if we have our lists defined
        if len(xerr) != len(x):
            xerr = [0] * len(x)

        if len(yerr) != len(y):
            yerr = [0] * len(y)

        if len(size) != len(x):
            size = [2] * len(x)

        if len(color) != len(x):
            color = ['black'] * len(x)

        #Check our axis limits
        if min(x) < self.x_min:
            self.x_min = min(x)

        if min(y) < self.y_min:
            self.y_min = min(y)

        if max(x) > self.x_max:
            self.x_max = max(x)

        if max(y) > self.y_max:
            self.y_max = max(y)

        #Reformat the axis labels to have the same precision
        fig, ax = plt.subplots()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))

        #Force ticks inward
        ax.tick_params(which = 'both', direction = 'in')
        
        #Create the scatterplot
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', markersize='0', ecolor=color)
        plt.scatter(x, y, s=size, c=color)

        #Set the axis limits
        plt.xlim(self.x_min, self.x_max)
        plt.ylim(self.y_min, self.y_max)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.draw()

        #Label the graph
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        #Plot the origin lines
        plt.axhline(0, color='black', linestyle='--')
        plt.axvline(0, color='black', linestyle='--')

        #If we are measuring regression then print it on the graph
        if self.bool_regression == "true":

            #Calculate the regression
            slope, intercept, r_value, p_value, std_err = linregress(x, y)

            #Add the legend
            textstr = r"$r^2=$"+str(round(r_value**2, 2))+"\nN: "+str(len(x))

            #Calculate top left
            text_x = self.x_min + 0.5
            text_y = self.y_max - 0.5
            plt.text(text_x, text_y, textstr, fontsize=20, 
                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))

            #Plot the line
            lin_fit = lambda x_val: slope*x_val + intercept

            xfit = (int(self.x_min) - 1, int(self.x_max) + 1)
            yfit = (lin_fit(int(self.x_min) - 1), lin_fit(int(self.x_max) + 1))  
            
            plt.plot(xfit, yfit, 'r', linestyle=':')

        #Show 1 to 1 line
        if self.onetoone == "true":
            plt.plot([-40, 40], [-40, 40], 'b', linestyle=':')

        #Check for SD Bounds
        if self.sd_bounds is not None:            
            plt.axhline(y=-self.sd_bounds, color='m', linestyle='--')
            plt.axhline(y=self.sd_bounds, color='m', linestyle='--')
            plt.axvline(x=-self.sd_bounds, color='m', linestyle='--')
            plt.axvline(x=self.sd_bounds, color='m', linestyle='--')

        #Set the size and dpi of the figure
        fig.set_size_inches(6, 6, forward=True)
        fig.set_dpi(200)

        #Need to not cut off half of our label
        plt.tight_layout()

        #Save the figure as a png
        plt.savefig(self.directory + self.output_prefix + "_scatter.png")
        plt.savefig(self.directory + self.output_prefix + "_scatter.ps")

        #Show if we are not running headless
        if self.pact_headless != "true":
            plt.show()

        return

    def xy_scatter(self, x, y, color = [], size = []):
        """Plot the dist of the read counts"""

        #Import the graphing toolkit
        from matplotlib import use, rcParams

        #Set the matplotlib base settings (headless, font)
        if self.pact_headless == "true":
            use('Agg')

        rcParams['font.sans-serif']='Arial'
        rcParams['font.size']='22'

        from matplotlib import pyplot as plt
        from matplotlib.ticker import FormatStrFormatter

        #Check if we have our lists defined
        if len(size) != len(x):
            size = [1] * len(x)

        if len(color) != len(x):
            color = ['black'] * len(x)

        #Check our axis limits
        if min(x) < self.x_min:
            self.x_min = min(x)

        if min(y) < self.y_min:
            self.y_min = min(y)

        if max(x) > self.x_max:
            self.x_max = max(x)

        if max(y) > self.y_max:
            self.y_max = max(y)

        #Reformat the axis labels to have the same precision
        fig, ax = plt.subplots()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))

        #Force ticks inward
        ax.tick_params(which = 'both', direction = 'in')
        
        #Create the scatterplot
        plt.scatter(x, y, s=size, color=color)

        #Set the axis limits
        plt.xlim(self.x_min, self.x_max)
        plt.ylim(self.y_min, self.y_max)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.draw()

        #Label the graph
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)

        #Plot the origin lines
        plt.axhline(0, color='black', linestyle='--')
        plt.axvline(0, color='black', linestyle='--')

        #If we are measuring regression then print it on the graph
        if self.bool_regression == "true":

            #Calculate the regression
            slope, intercept, r_value, p_value, std_err = linregress(x, y)

            #Add the legend
            textstr = r"$r^2=$"+str(round(r_value**2, 2))+"\nN: "+str(len(x))

            #Calculate top left
            text_x = self.x_min + 0.5
            text_y = self.y_max - 0.5
            plt.text(text_x, text_y, textstr, fontsize=20, 
                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white'))

            #Plot the line
            lin_fit = lambda x_val: slope*x_val + intercept

            xfit = (int(self.x_min) - 1, int(self.x_max) + 1)
            yfit = (lin_fit(int(self.x_min) - 1), lin_fit(int(self.x_max) + 1))  
            
            plt.plot(xfit, yfit, 'r', linestyle=':')

        #Show 1 to 1 line
        if self.onetoone == "true":
            plt.plot([-40, 40], [-40, 40], 'b', linestyle=':')

        #Check for SD Bounds
        if self.sd_bounds is not None:            
            plt.axhline(y=-self.sd_bounds, color='m', linestyle='--')
            plt.axhline(y=self.sd_bounds, color='m', linestyle='--')
            plt.axvline(x=-self.sd_bounds, color='m', linestyle='--')
            plt.axvline(x=self.sd_bounds, color='m', linestyle='--')

        #Set the size and dpi of the figure
        fig.set_size_inches(6, 6, forward=True)
        fig.set_dpi(200)

        #Need to not cut off half of our label
        plt.tight_layout()

        #Save the figure as a png
        plt.savefig(self.directory + self.output_prefix + "_scatter.png")
        plt.savefig(self.directory + self.output_prefix + "_scatter.ps")

        #Show if we are not running headless
        if self.pact_headless != "true":
            plt.show()

        return

    def create_csv(self, x, y):
        """Output a CSV for other programs"""

        #Open the output file and write it
        with open(self.directory + self.output_prefix + "_x_vs_y.csv", "w") as csv_out:

            #Write the header
            csv_out.write(self.x_label + "," + self.y_label + "\n")

            #Create the output lines
            csv_out.write('\n'.join([str(xlist) + "," + str(ylist) for xlist,ylist in zip(x, y)]))
        return

    def shared_mutations(self, dict_merged_datasets):
        """Return a list of amino acid points"""

        countx = 0
        #Count the basal winners in set x
        for key in dict_merged_datasets[self.dataset_x]:
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN":
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Count our basal winner
                if dict_merged_datasets[self.dataset_x][key][mut][self.x_column] >= self.winner_threshold:
                    countx += 1

        county = 0
        #Count the basal winners in set y
        for key in dict_merged_datasets[self.dataset_y]:
            for mut in dict_merged_datasets[self.dataset_y][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN":
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Count our basal winner
                if dict_merged_datasets[self.dataset_y][key][mut][self.y_column] >= self.winner_threshold:
                    county += 1

        #Count the shared mutations in x and y
        countxy = 0
        #Loop the locations
        for key in dict_merged_datasets[self.dataset_x]:

            #Check if the location is not in y, skip if so
            if key not in dict_merged_datasets[self.dataset_y]:
                continue

            #Loop the muts
            for mut in dict_merged_datasets[self.dataset_x][key]:
                
                #Skip WT
                if mut == dict_merged_datasets[self.dataset_x][key][mut]['wt_residue']:
                    continue

                #Skip WT
                if mut == dict_merged_datasets[self.dataset_y][key][mut]['wt_residue']:
                    continue

                #Skip NaN
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] == "NaN" 
                    or dict_merged_datasets[self.dataset_y][key][mut][self.y_column] == "NaN"):
                    continue

                #See if we are above our threshold
                if (dict_merged_datasets[self.dataset_x][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_ref_counts"] < self.ref_threshold or
                    dict_merged_datasets[self.dataset_x][key][mut]["adjusted_sel_counts"] < self.sel_threshold or
                    dict_merged_datasets[self.dataset_y][key][mut]["adjusted_sel_counts"] < self.sel_threshold):
                    continue

                #Count our basal winner
                if (dict_merged_datasets[self.dataset_x][key][mut][self.x_column] >= self.winner_threshold and
                    dict_merged_datasets[self.dataset_y][key][mut][self.y_column] >= self.winner_threshold):
                    countxy += 1

        return countx, county, countxy

    def set_vs_set(self, dict_merged_datasets, dict_color={}, dict_size={}, dict_error={}):
        """Main entrypoint"""

        #Check to see if our dataset is in our merged datasets
        if (self.dataset_x not in dict_merged_datasets or
            self.dataset_y not in dict_merged_datasets):
            print("[Set vs Set Error] Check the config file, the dataset wanted was not combined.")
            quit()

        #Check if we have a custom point size
        if len(dict_size) > 0:
            #Parse the size dict to list
            size = self.dictsize_to_list(dict_merged_datasets, dict_size)
        else:
            size = []

        #Check if we have a custom errorbars
        if len(dict_error) > 0:
            #Parse the error dict to list
            xerr, yerr = self.dicterror_to_list(dict_merged_datasets, dict_error)
        else:
            xerr = []
            yerr = []

        #Select which color table to use
        str_color_table = "" #String to return to

        if self.point_color == "outlier":
            #Find outliers
            color, str_color_table = self.outliers(dict_merged_datasets)            
        elif self.point_color == "winner":
            #Find winners
            color, str_color_table = self.winners(dict_merged_datasets)
        elif self.point_color == "outlier_sign":
            #Find outliers color by sign diff
            color, str_color_table = self.outlier_sign(dict_merged_datasets)
        elif self.point_color == "amino":
            #Find aminos
            color = self.aminos(dict_merged_datasets)           
        elif self.point_color == "classifier_color":
            #Use an outside dict
            color = self.dictcolor_to_list(dict_merged_datasets, dict_color)
        else:
            color = []

        #Get the two lists of values
        if self.xy_grouping == "grouped_location_outlier":
            #Group by location
            x, y, x_nan, y_nan, xerr, yerr, color, size = self.xy_list_grouplocout(dict_merged_datasets, dict_color)
        else:
            #The standard point
            x, y, x_nan, y_nan = self.xy_list(dict_merged_datasets)
            
        #Select which XY Scatter Plotter To Use
        if self.xy_type == "standard":
            #Call the standard xy plotter
            self.xy_scatter(x, y, color, size)
        elif self.xy_type == "errorbars":
            #Call the error bar scatter plotter
            self.xy_scatter_error(x, y, color, size, xerr, yerr)

        #Output a CSV of the data
        if self.output_csv == "true":
            self.create_csv(x, y)

        #Perform regression?
        if self.bool_regression == "true":
            #Perform the regression
            slope, intercept, r_value, p_value, std_err = linregress(x, y)

            #Print regression stats
            str_regression = '\n'.join([
                "Regression for the entire dataset:",
                "slope = " + str(round(slope,3)),
                "intercept = " + str(round(intercept,3)),
                "r = " + str(round(r_value,3)),
                "r^2 = " + str(round(r_value**2,3)),
                "N = " + str(len(x)),
                ])
            print(str_regression)
        else:
            str_regression = ""

        #Do we calculate shared winners counts
        if self.count_shared == "true":
            #Get our basal counts of shared winners
            countx, county, countxy = self.shared_mutations(dict_merged_datasets)

            #Print shared winner stats
            str_shared = '\n'.join([
                "Dataset vs Dataset",
                "Dataset X: " + self.dataset_x,
                "Dataset Y: " + self.dataset_y,
                "Counts of nonsynonymous mutations above threshold: " + str(self.winner_threshold),
                "Dataset X: " + str(countx),
                "Dataset Y: " + str(county),
                "Shared: " + str(countxy),
                ])
            print(str_shared)
        else:
            str_shared = ""

        return str_color_table + "\n" + str_regression  + "\n" + str_shared

if __name__ == '__main__':
    #Remind the user that the protocol needs to be ran within the context of PACT
    print("[Set vs Set Error] This file needs to be ran within the context of PACT.")
