#!/usr/bin/env python
# -*-coding: utf-8 -*-

"""
Usage: $ {1: program}.py
"""
import json
import re

import db

db_objects = db.get_db()
reactome_db = db_objects["reactome_db"]
reactome_db_name = db_objects["reactome_db_name"]
reactome_coll = db_objects["reactome_coll"]
reactome_coll_name = db_objects["reactome_coll_name"]
reactome_aql = db_objects["reactome_aql"]


def clean_label(label):
    """Add quotes or escape a double quote"""
    label.replace('"', '"')
    if " " in label:
        label = f'"{label}"'
    elif ")" in label:
        label = f'"{label}"'

    return label


def ts_convert(ts):
    """Convert timestamp from Reactome format to ISO format"""

    ts = re.sub(r"(\d{4,4}-\d{2,2}-\d{2,2})\s(\d{2,2}:\d{2,2}:\d{2,2}).*", r"\1T\2.000Z", ts)
    return ts


def stats():

    query = f"""
    FOR d IN {reactome_coll_name}
        FILTER d.doc_type == "Reaction"
        RETURN d
    """

    s = {"class_names": {}, "schema_classes": {}, "categories": {}}

    cursor = reactome_aql.execute(query)
    for c in cursor:
        class_name = c["doc"]["className"]
        schema_class = c["doc"]["schemaClass"]
        category = c["doc"].get("category")

        if class_name in s["class_names"]:
            s["class_names"][class_name] += 1
        else:
            s["class_names"][class_name] = 0

        if schema_class in s["schema_classes"]:
            s["schema_classes"][schema_class] += 1
        else:
            s["schema_classes"][schema_class] = 0

        if category in s["categories"]:
            s["categories"][category] += 1
        else:
            s["categories"][category] = 0

    import json

    print("Stats:\n", json.dumps(s, indent=4))


def main():
    stats()


if __name__ == "__main__":
    main()
