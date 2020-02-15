import datetime
import json
import os.path
import time

import requests
import structlog

import db
import settings
from settings import REACTOME_API_URL, appdir

log = structlog.getLogger(__name__)

db_objects = db.get_db()
reactome_db = db_objects["reactome_db"]
reactome_db_name = db_objects["reactome_db_name"]
reactome_coll = db_objects["reactome_coll"]
reactome_coll_name = db_objects["reactome_coll_name"]

start_time = datetime.datetime.now()
downloaded_count = 0


# TODO - rework this entire file
#    1. Given species - get_starting_events
#    2. Recursively process each event (pathways have children, reactions are nodes)
#    3. Check if already in db, if not save event in arangodb
#    Collect additional (dbid) records, e.g. proteins, etc?


def get_starting_events(species_id):
    """Get top event (generally top pathways) for given species"""

    url = f"{reactome_api_url}/eventsHierarchy/{species_id}"
    r = requests.get(url)

    if r.status_code != 200:
        log.error(f"Could not get pathway hierarchy for {species_id}")
        quit()

    events = r.json()

    return events


def event_in_database(stid):
    """Check to see if stId is in database"""

    query = f"""
    FOR d in {reactome_coll_name}
        d.doc.stID == "{stid}"
        RETURN d
    """
    cursor = reactome_db.aql.execute(query)
    if len(list(cursor)) > 0:
        return True

    return False


def save_event_to_database(event):
    """Save to database"""

    # TODO check format of events in database and make sure this matches
    pass


def get_event_from_database(stid):
    """Get event from database"""

    pass


def get_event(stid: str):
    """Get pathway or reaction data structure and save it"""

    if not event_in_database(stid):
        time.sleep(0.2)
        r = requests.get(f"{reactome_api_url}/query/enhanced/{stid}")
        if r.status_code != 200:
            log.error("Could not retrieve {stid}")
            return True

        event = r.json()

        save_event_to_database(event)

    else:
        event = get_event_from_database(stid)

    save
    return event


def process_children(children: list):
    """Recursively process hierarchical response child elements"""

    for event in children:
        if "stId" in event:
            get_event(event["stId"])
        if "children" in event:
            process_children(event["children"])
