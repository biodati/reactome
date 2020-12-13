# Standard Library
import inspect
import logging
import os
import os.path

# Third Party Imports
import bel.core.settings as bel_settings

# Reactome to BEL Settings

REACTOME_API_URL = "https://reactome.org/ContentService/data"

# https://reactome.org/ContentService/data/query/enhanced/390750

# app directory
appdir = os.path.split(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))[0]
ARANGO_DATA_DIR = f"{appdir}/data/arangodb"

os.makedirs(ARANGO_DATA_DIR, exist_ok=True)


ARANGODB = "http://thor:9529"
REACTOME_DB = "_system"
REACTOME_COLL = "reactome_2020q4"
# REACTOME_COLL = "reactome"

LOGLEVEL = logging.INFO


# BEL Settings
bel_settings.ARANGO_URL = "http://dev.biodati.test:8529"
bel_settings.ELASTICSEARCH_URL = "http://dev.biodati.test:9200"


amino_acids = {
    "alanine": "ala",
    "arginine": "arg",
    "asparagine": "asn",
    "aspartic acid": "asp",
    "cysteine": "cys",
    "glutamic acid": "glu",
    "glutamine": "gln",
    "glycine": "gly",
    "histidine": "his",
    "isoleucine": "ile",
    "leucine": "leu",
    "lysine": "lys",
    "methionine": "met",
    "phenylalanine": "phe",
    "proline": "pro",
    "serine": "ser",
    "threonine": "thr",
    "tryptophan": "trp",
    "tyrosine": "tyr",
    "valine": "val",
    "unknown": "?",
}

amino_acids_regex = "|".join([aa for aa in sorted(amino_acids.keys(), key=len)])
