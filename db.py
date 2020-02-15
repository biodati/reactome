# Third Party Imports
import time
from typing import List

import arango
import requests
import structlog

# Local Imports
import settings

log = structlog.getLogger(__name__)

reactome_db_name = settings.REACTOME_DB
reactome_coll_name = settings.REACTOME_COLL

try:
    client = arango.ArangoClient(hosts=settings.ARANGODB)
except Exception as e:
    log.error(f"Could not get client handle to arangodb: {e}")
    quit()


def get_db() -> arango.database.StandardDatabase:
    """Get arangodb database handle

    Args:

    Returns:
        arango.database.StandardDatabase: Description
    """

    # client is created when module is first imported
    sys_db = client.db("_system")

    # Create a new database named "nanopubs"
    if sys_db.has_database(reactome_db_name):
        reactome_db = client.db(reactome_db_name)
    else:
        sys_db.create_database(name=reactome_db_name)
        reactome_db = client.db(reactome_db_name)

    # Add nanopubs collection
    if reactome_db.has_collection(reactome_coll_name):
        reactome_coll = reactome_db.collection(reactome_coll_name)
    else:
        reactome_coll = reactome_db.create_collection(reactome_coll_name, index_bucket_count=64)
        # Adding of hash indexes:
        reactome_coll.add_hash_index(fields=["doc_type"], name="doc_type_idx")
        reactome_coll.add_hash_index(fields=["doc.stId"], name="stid_idx")

    reactome_aql = reactome_db.aql
    reactome_aql.cache.configure(mode="on", max_results=10000)

    return {
        "reactome_db": reactome_db,
        "reactome_db_name": reactome_db_name,
        "reactome_coll": reactome_coll,
        "reactome_coll_name": reactome_coll_name,
        "reactome_aql": reactome_aql,
    }


# #############################################################################
# Initialize arango_client !!!!!!!!!!!!!!!!!!!
# #############################################################################

db_objects = get_db()
reactome_db = db_objects["reactome_db"]
reactome_db_name = db_objects["reactome_db_name"]
reactome_coll = db_objects["reactome_coll"]
reactome_coll_name = db_objects["reactome_coll_name"]
reactome_aql = db_objects["reactome_aql"]


def reset_database():
    """Reset the nanopub database"""

    # client is created when module is first imported
    sys_db = client.db("_system")
    try:
        sys_db.delete_database(settings.REACTOME_DB)
    except Exception as e:
        pass

    get_db()


def save_entity(doc):

    doc_type = doc["className"]
    rid = doc["dbId"]

    reactome_coll.insert({"_key": str(rid), "doc_type": doc_type, "doc": doc})


def has_stid(stid):
    """Check if stid exists in arangodb"""

    query = f"""
    FOR d in {reactome_coll_name}
        FILTER d.doc.stId == "{stid}"
        RETURN d._key
    """

    cursor = reactome_aql.execute(query)
    return list(cursor)


def get_entity(dbid):
    """Get entity and save it"""

    dbid = str(dbid)

    # print('dbId', dbid)

    if reactome_coll.has(dbid):
        return reactome_coll.get(dbid)["doc"]

    else:
        time.sleep(0.1)
        r = requests.get(f"{settings.REACTOME_API_URL}/query/enhanced/{dbid}")
        if r.status_code != 200:
            log.error(f"Could not retrieve {dbid}")
            return True

        doc = r.json()
        save_entity(doc)
        return doc


def initial_reload():
    """Used to reload downloaded content into arangodb"""

    import glob
    import json

    reset_database()

    files = glob.glob(f"./data/9606/*.json")
    for fn in files:
        with open(fn, "r") as f:
            doc = json.load(f)

        doc_type = doc["className"]
        species = "9606"
        rid = doc["dbId"]

        # print(f"_key: {rid}   doc_type: {doc_type}   doc: {doc['dbId']}")

        reactome_coll.insert({"_key": str(rid), "doc_type": doc_type, "doc": doc})


def collect_all_reactions() -> List[str]:
    """Get all reaction dbId's"""

    query = f"""
    FOR d in {reactome_coll_name}
        FILTER d.doc_type == "Reaction"
        RETURN d._key
    """
    cursor = reactome_aql.execute(query, batch_size=10, ttl=3600000)

    return list(cursor)
