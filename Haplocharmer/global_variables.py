from pathlib import Path

DOCS = "https://Haplocharmer.readthedocs.io/en/latest/"
GIT_URL = "https://github.com/thdurand4/HaploCharmer"

INSTALL_PATH = Path(__file__).resolve().parent
#SINGULARITY_URL_FILES = [('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.culebront_tools.sif',
#                          'INSTALL_PATH/containers/Singularity.culebront_tools.sif'),
#                         ('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.report.sif',
#                          'INSTALL_PATH/containers/Singularity.report.sif')
#                         ]

DATATEST_URL_FILES = ("https://itrop.ird.fr/culebront_utilities/Data-Xoo-sub.zip", "Data-Xoo-sub.zip")
SCRIPTS = INSTALL_PATH.joinpath("scripts")



