#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path

from snakemake.logging import logger
from snakemake.utils import validate
import re
import os.path
from .global_variables import *
from snakecdysis import *

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



class HaploCharmer(SnakEcdysis):
    """
    to read file config
    """

    def __init__(self, dico_tool, workflow, config):
        super().__init__(**dico_tool, workflow=workflow, config=config)
        self.config = config
        self.use_env_modules = workflow.use_env_modules
        self.use_singularity = workflow.use_singularity
        # workflow is available only in __init__
        # print("\n".join(list(workflow.__dict__.keys())))
        # print(workflow.__dict__)

        # Initialisation of Haplocharmer attributes
        
        self.__check_config_dic()
        
        #self.write_config(self.config)
        
    def __check_config_dic(self):
        """Configuration file checking"""
        self.__get_tools_config("ENVMODULE")
       
        #self.tools_activated = self.__build_tools_activated("OPTIONAL", ("REMOVE_DUPLICATES"), True)
        
        self._check_dir_or_string(level1="DATA", level2="OUTPUT")
        self._check_dir_or_string(level1="DATA", level2="SCRIPTS")
        
        if (os.path.isfile(self.get_config_value(level1='DATA', level2='REFERENCE'))):
            pass
        else:
            raise ValueError(f"{bcolors.BOLD}{bcolors.FAIL}File {bcolors.OKBLUE}{self.get_config_value(level1='DATA', level2='REFERENCE')} {bcolors.FAIL}doesn't exist")
        
        if (os.path.isfile(self.get_config_value(level1='DATA', level2='BEDFILE'))):
            pass
        else:
            raise ValueError(f"{bcolors.BOLD}{bcolors.FAIL}File {bcolors.OKBLUE}{self.get_config_value(level1='DATA', level2='BEDFILE')} {bcolors.FAIL}doesn't exist")
        
        if (os.path.isfile(self.get_config_value(level1='DATA', level2='INFORMATION_FILE'))):
            pass
        else:
            raise ValueError(f"{bcolors.BOLD}{bcolors.FAIL}File {bcolors.OKBLUE}{self.get_config_value(level1='DATA', level2='INFORMATION_FILE')} {bcolors.FAIL}doesn't exist")
        if (self.get_config_value(level1='OPTIONAL', level2='REMOVE_DUPLICATES') == True or self.get_config_value(level1='OPTIONAL', level2='REMOVE_DUPLICATES') == False):
            pass
        else:
            raise ValueError(f"{bcolors.BOLD}{bcolors.FAIL}ERROR {bcolors.OKBLUE}{self.get_config_value(level1='OPTIONAL', level2='REMOVE_DUPLICATES')} {bcolors.FAIL}is not allowed. Please make sure to write : {bcolors.OKGREEN}True {bcolors.FAIL}or {bcolors.OKGREEN}False {bcolors.FAIL}in {bcolors.OKGREEN}REMOVE_DUPLICATES {bcolors.FAIL}section")
        
        def chek_info_file(self, path):
            # check of header
            infile = open(path, 'r')
            if not 'R1\tR2\tSM\tLB\tID\tPU\tPLATFORM' in infile.readline():
                raise ValueError(
                    f"{bcolors.BOLD}{bcolors.FAIL}INFORMATION FILE CHECKING FAIL : Please add {bcolors.OKGREEN}R1\tR2\tSM\tLB\tID\tPU\tPLATFORM {bcolors.BOLD}{bcolors.FAIL}in INFORMATION FILE tabulated header !")
            infile.close()
        
        self.info = chek_info_file(self, self.get_config_value(level1='DATA', level2='INFORMATION_FILE'))
        
        


    def __check_tools_config(self, tool, mandatory=[]):
        """Check if path is a file if not empty
        :return absolute path file"""
    
        # If only envmodule
        if self.use_env_modules and not self.use_singularity:
            envmodule_key = self.tools_config["ENVMODULE"][tool]
            if not envmodule_key:
                raise ValueError(
                    f'{bcolors.BOLD}{bcolors.FAIL}CONFIG FILE CHECKING FAIL : please check {bcolors.OKGREEN}tools_config.yaml {bcolors.BOLD}{bcolors.FAIL}in the {bcolors.OKGREEN}"ENVMODULE" {bcolors.BOLD}{bcolors.FAIL}section, {bcolors.OKBLUE}{tool} {bcolors.BOLD}{bcolors.FAIL}is empty')
        

        # If envmodule and singularity
        if self.use_env_modules and self.use_singularity:
            raise ValueError(
                f"{bcolors.BOLD}{bcolors.FAIL}CONFIG FILE CHECKING FAIL : Use env-module or singularity but don't mix them")

    
    def __get_tools_config(self, key):
    	for tool, activated in self.tools_config[key].items():
    	    if tool != "None":
       	    	self.__check_tools_config(tool, [tool])
    
    
