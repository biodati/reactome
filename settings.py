import inspect
import os
import os.path

REACTOME_API_URL = "https://reactome.org/ContentService/data"

# app directory
appdir = os.path.split(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))[0]
ARANGO_DATA_DIR = f"{appdir}/data/arangodb"
try:
    os.makedirs(ARANGO_DATA_DIR)
except Exception as e:
    pass

ARANGODB = "http://thor:9529"
REACTOME_DB = "_system"
REACTOME_COLL = "reactome"
