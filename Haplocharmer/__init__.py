#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Haplocharmer.global_variables import *
from pathlib import Path
from .global_variables import GIT_URL, DOCS, DATATEST_URL_FILES, SINGULARITY_URL_FILES

logo = Path(__file__).parent.resolve().joinpath('Haplocharmer_logo.png').as_posix()

__version__ = Path(__file__).parent.resolve().joinpath("VERSION").open("r").readline().strip()


__doc__ = """HaploCharmer   is   a   Snakemake   pipeline   for   small   haplotypes   calling   from   short
read   sequencing   data.   This   pipeline   is   aimed   to   provide   more   informative
genotyping   than   the   standard   SNPs,   particularly   for   polyploid   species,   in   the
context   of   quantitative   genetics,   population   genetics   and   genetic   mapping
analyses."""

description_tools = f"""
    Welcome to Haplocharmer version: {__version__} ! Created on May 2023
    @author: Théo Durand (CIRAD)
    @email: theo.durand@cirad.fr

    Please cite our github: GIT_URL
    and GPLv3 Intellectual property belongs to CIRAD and authors.
    Documentation avail at: DOCS"""

dico_tool = {
    "soft_path": Path(__file__).resolve().parent.as_posix(),
    "url": GIT_URL,
    "docs": DOCS,
    "description_tool": description_tools,
    "singularity_url_files": SINGULARITY_URL_FILES,
    "datatest_url_files": DATATEST_URL_FILES,
    "snakefile": Path(__file__).resolve().parent.joinpath("snakefiles", "Snakefile"),
    "snakemake_scripts": Path(__file__).resolve().parent.joinpath("snakemake_scripts")
}
