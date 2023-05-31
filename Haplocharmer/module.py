#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path

from snakemake.logging import logger
from snakemake.utils import validate
import re
from .global_variables import *
from snakecdysis import *


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
        
        self.write_config(f"{self.config['DATA']['OUTPUT']}/config_corrected.yaml")
        
    def __check_config_dic(self):
        """Configuration file checking"""

        self.tools_activated = self.__build_tools_activated("OPTIONAL", ("REMOVE_DUPLICATES"), True)
        self._check_dir_or_string(level1="DATA", level2="REFERENCE")
        self._check_dir_or_string(level1="DATA", level2="OUTPUT")
        self._check_dir_or_string(level1="DATA", level2="SCRIPTS")
        self._check_dir_or_string(level1="DATA", level2="BEDFILE")
        self._check_dir_or_string(level1="DATA", level2="INFORMATION_FILE")

    def __check_tools_config(self, tool, mandatory=[]):
        """Check if path is a file if not empty
        :return absolute path file"""
        tool_OK = True

        # If only envmodule
        if self.use_env_modules and not self.use_singularity:
            envmodule_key = self.tools_config["ENVMODULE"][tool]
            if not envmodule_key:
                raise ValueError(
                    f'CONFIG FILE CHECKING FAIL : please check tools_config.yaml in the "ENVMODULE" section, {tool} is empty')
            tool_OK = True

        # If envmodule and singularity
        if self.use_env_modules and self.use_singularity:
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL : Use env-module or singularity but don't mix them")

        if len(mandatory) > 0 and not tool_OK:
            raise FileNotFoundError(
                f'CONFIG FILE CHECKING FAIL : please check tools_config.yaml in the  {tool} params, please append Singularity or module load, is mandatory for tool: {" ".join(mandatory)}')
    
    
    
    def __build_tools_activated(self, key, allow, mandatory=False):
        tools_activate = []
        for tool, activated in self.config[key].items():
            if tool in allow:
                boolean_activated = var_2_bool(key, tool, activated)
                if boolean_activated:
                    tools_activate.append(tool)
                    self.config[key][tool] = boolean_activated
                    self.__check_tools_config(tool, [tool])
            else:
                raise ValueError(f'CONFIG FILE CHECKING FAIL : {key} {tool} not allowed on HaploCharmer"')
        if len(tools_activate) == 0 and mandatory:
            raise ValueError(f"CONFIG FILE CHECKING FAIL : you need to set True for at least one {key} from {allow}")
        return tools_activate


    def __var_2_bool(self, key, tool, to_convert):
        """convert to boolean"""
        if isinstance(type(to_convert), bool):
            return to_convert
        elif f"{to_convert}".lower() in ("yes", "true", "t"):
            return True
        elif f"{to_convert}".lower() in ("no", "false", "f"):
            return False
        else:
            raise TypeError(
                f'CONFIG FILE CHECKING FAIL : in the "{key}" section, "{tool}" key: "{to_convert}" is not a valid boolean')

