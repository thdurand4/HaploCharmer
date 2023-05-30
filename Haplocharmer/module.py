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
        # workflow is available only in __init__
        # print("\n".join(list(workflow.__dict__.keys())))
        # print(workflow.__dict__)

        # Initialisation of Haplocharmer attributes
