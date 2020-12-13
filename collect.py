#!/usr/bin/env python
# -*-coding: utf-8 -*-

"""
Usage: $ {1: program}.py
"""

# Standard Library
import datetime
import json
import os.path
import time

# Third Party Imports
import requests

# Local Imports
import db
import settings
from logging_setup import log
from settings import REACTOME_API_URL, appdir

reactome_db = db.reactome_db
reactome_db_name = db.reactome_db_name
reactome_coll = db.reactome_coll
reactome_coll_name = db.reactome_coll_name

start_time = datetime.datetime.now()
downloaded_count = 0

"""Overview

This script will download the species (human=9606) pathways from Reactome
and all of the reactions part of those pathways recursively in the 
sub-pathways.

The main.py script will then pull any entities referenced in the reactions that are
needed to process Reactome into BEL.
"""


def get_starting_entities(species_id):
    """Get top entities (generally top pathways) for given species"""

    url = f"{REACTOME_API_URL}/eventsHierarchy/{species_id}"
    r = requests.get(url)

    if r.status_code != 200:
        log.error(f"Could not get pathway hierarchy for {species_id}")
        quit()

    events = r.json()

    return events


def save_entity(stid: str):
    """Get pathway or reaction data structure and save it"""

    global downloaded_count

    if not db.is_stid_in_database(stid):
        downloaded_count += 1

        print(f"Downloading {stid}  Cnt: {downloaded_count}")

        time.sleep(0.2)

        url = f"{REACTOME_API_URL}/query/enhanced/{stid}"
        r = requests.get(url)
        if r.status_code == 200:
            event = r.json()
            db.save_entity(event)
        else:
            print(f"Could not retrieve {stid} at {url}")


def process_entities(entities: list):
    """Recursively process hierarchical response child elements"""

    for entity in entities:
        if "stId" in entity:
            save_entity(entity["stId"])
            print("STID", entity["stId"])

        if "children" in entity:
            process_entities(entity["children"])


def main():

    entities = get_starting_entities(9606)
    process_entities(entities)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted")
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
