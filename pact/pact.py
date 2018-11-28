#!/usr/bin/env python3
#
# Protein Analysis and Classifier Toolkit
# Author: Justin R. Klesmith
# Copyright (C) 2018 by Regents of the University of Minnesota
# Copyright (C) 2018 by Justin R. Klesmith
#
# This software is released under GNU General Public License 3
# Additional license options available at http://license.umn.edu
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

"""Entry point wrapper for the Protein Analysis and Classifier Toolkit"""

from sys import version_info

#Check to see if our version is supported
req_version = (3,4)
cur_version = version_info

if cur_version < req_version:
    print("[PACT Error] Your Python interpreter is too old."
          "Minimium version required is Python 3.4 (64 bit recommended).")
    quit()

from os.path import isfile, dirname, realpath
from sys import platform, modules, argv
from time import strftime
from configparser import ConfigParser, NoSectionError, NoOptionError
from importlib import import_module
from multiprocessing import freeze_support
from pact.pact_common import file_checker, command_line_args

#Set the author information
__author__ = "Justin R. Klesmith and Benjamin J. Hackel"
__copyright__ = ["Copyright (C) 2018 by Regents of the University of Minnesota", 
                 "Copyright (C) 2018 Justin R. Klesmith"]
__license__ = "GPL-3.0"
__version__ = "2018.6"
__maintainer__ = "Justin R. Klesmith"
__email__ = ["jrk@umn.edu", "hackel@umn.edu"]

#Get the PACT script path
pact_path = dirname(realpath(argv[0]))

#Check to see if the protocols links exist
if file_checker(pact_path + "/pact/pact_protocols.ini"):
    #Load the module links config file into a directory
    config_parser = ConfigParser()

    #Load the global and workflow config elements
    config_parser.read(pact_path + "/pact/pact_protocols.ini")

    #Loop the module
    try:
        dict_protocols = {mapping[0].lower(): mapping[1].lower() 
                          for mapping in config_parser.items("protocols")}
    except NoSectionError:
        print("[PACT Error] The config file pact_protocol.ini is incorrect.")
        print("[PACT Error] There is something wrong with the [protocols] naming of the file.")
        quit()
    except NoOptionError:
        print("[PACT Error] The config file pact_protocol.ini is incorrect.")
        print("[PACT Error] There is something wrong with the name of a option flag.")
        quit()

#Check to see if the program links exist
if file_checker(pact_path + "/pact/pact_external_programs.ini"):
    #Load the global and workflow config elements
    config_parser.read(pact_path + "/pact/pact_external_programs.ini")

    #Loop the module
    try:   
        dict_programs = {mapping[0].lower(): mapping[1] 
                         for mapping in config_parser.items("programs")}
    except NoSectionError:
        print("[PACT Error] The config file pact_external_programs.ini is incorrect.")
        print("[PACT Error] There is something wrong with the [programs] naming of the file.")
        quit()
    except NoOptionError:
        print("[PACT Error] The config file pact_external_programs.ini is incorrect.")
        print("[PACT Error] There is something wrong with the name of a option flag.")
        quit()

    #Check to see if each file exists
    for mapping in dict_programs:
        if not isfile(dict_programs[mapping]):
            print("[PACT Error] The system cannot find the program "+mapping+" at "+dict_programs[mapping])

#Parse the command line for inputs
file_config = command_line_args('-c')
if file_config[0]:
    file_checker(file_config[1])
else:
    print("[PACT Error] The command line for the config file -c was not found.")
    quit()

#Wrap all of the code into main
if __name__ == '__main__':
    """Let's wrap all of our imported code in the main module."""
    #If we are on windows we need to have freeze support for multithreading
    if platform.startswith('win'):
        #Not required if called from __main__ but the modules may be imported from a different script
        if __name__ != "__main__":
            modules["__main__"] = modules[__name__]
        freeze_support()

    #Open our config file and see which protocol we want
    with open(file_config[1]) as input_config:
        for line in input_config.readlines():
            #Get the pact config version
            if line.startswith("pact_config_version:"):
                config_version = line.split(":")[1].lower().rstrip("\n\r ").lstrip(" ")

            #Get the protocol
            if line.startswith("pact_protocol:"):
                protocol = line.split(":")[1].lower().rstrip("\n\r ").lstrip(" ")

    #Print the preamble    
    string_preamble = ("[PACT] Protein Analysis and Classifier Toolkit\n"
                       "[PACT] "+__copyright__[0]+"\n"
                       "[PACT] "+__copyright__[1]+"\n"
                       "[PACT] Author: "+__author__+"\n"
                       "[PACT] Contact: "+__email__[0]+", "+__email__[1]+"\n"
                       "[PACT] Version: "+__version__+"\n"
                       "[PACT] License: "+__license__+"\n"
                       "[PACT] Please cite:\n"
                       "[PACT] Github [user: JKlesmith] (https://github.com/JKlesmith/PACT)\n"
                       "[PACT] \n"
                       "[PACT] "+strftime("%H:%M:%S")+"\n"
                       "[PACT] "+strftime("%m/%d/%Y")+"\n"
                       "[PACT] PACT Path: " + pact_path + "\n"
                       "[PACT] Config File: "+file_config[1]+"\n"
                       "[PACT] \n"       
                       "[PACT] Loading Protocol: " + protocol)
    print(string_preamble)

    #Check to see if we can import our module
    if protocol not in dict_protocols:
        print("[PACT Error] Cannot find the requested protocol.")
        quit()

    #Load our protocol
    try:
        loaded_protocol = import_module(dict_protocols[protocol])
    except ImportError:
        print("[PACT Error] Cannot find the requested protocol.")
        quit()

    #Get the major class from our protocol
    try:
        protocol_class = getattr(loaded_protocol, protocol)
    except AttributeError:
        print("[PACT Error] Cannot find the main class for the protocol.")
        quit()

    #Create our protocol object
    protocol = protocol_class(file_config[1], dict_programs, string_preamble)

    #Check our protocol version
    print("[PACT] Protocol version: " + protocol.version())
    if protocol.version() != __version__:
        print("[PACT Error] The version of the protocol does not match the PACT entrypoint.")
        quit()

    #Check our config version
    print("[PACT] Config file version: " + config_version)
    if config_version != __version__:
        print("[PACT Error] The version of the config file does not match the PACT entrypoint.")
        quit()

    #Load our main entrypoint
    protocol.protocol()

    #Complete
    print("[PACT] Run Complete")
