#!/usr/bin/env python
# -*-coding: utf-8 -*-

"""
Usage: $ {1: program}.py
"""

# Standard Library
import json
import os
import sys
import traceback

# Third Party Imports
import structlog
from logging_tree import printout

# Local Imports
import collect
import db
import settings
import to_bel
from logging_setup import log

db_objects = db.get_db()
reactome_coll_name = db_objects["reactome_coll_name"]
reactome_aql = db_objects["reactome_aql"]


def pathways_to_reactions(pathways: list = [], reactions: set = set()):
    """Get pathway reactions - only used if need to filter reactions processed to some pathways"""

    for pathway in pathways:
        filter = f'FILTER d.doc.stId == "{pathway}"'
        query = f"""
        FOR d in {reactome_coll_name}
            {filter}
            RETURN d
        """
        cursor = reactome_aql.execute(query, batch_size=10, ttl=3600)

        for doc in cursor:
            for reaction in doc["doc"]["hasEvent"]:
                # print("Reaction", reaction)
                if not isinstance(reaction, dict):
                    continue
                if reaction["className"] == "Reaction":
                    reactions.add(reaction["dbId"])
                elif reaction["className"] == "Pathway":
                    reactions = pathways_to_reactions(
                        pathways=[reaction["stId"]], reactions=reactions
                    )
                    # print(f"Pathway: {pathway}  ClassName: {reaction['className']}")

    return reactions


def review_reactions():
    reaction_ids = to_bel.collect_all_reactions()
    stats = {}
    for rid in reaction_ids:
        doc = db.get_entity(rid)
        dbid = doc["dbId"]
        schema_class = doc["schemaClass"]
        category = doc["category"]
        input_cnt = len(doc.get("input", []))
        output_cnt = len(doc.get("output", []))
        input_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("input", [])
            if isinstance(t, dict)
        ]
        output_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("output", [])
            if isinstance(t, dict)
        ]
        catalyst_cnt = len(doc.get("catalystActivity", []))
        catalyst_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("catalystActivity", [])
            if isinstance(t, dict)
        ]
        regulatedby_cnt = len(doc.get("regulatedBy", []))
        regulatedby_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("regulatedBy", [])
            if isinstance(t, dict)
        ]
        print(
            f"{dbid} class: {schema_class} category: {category} inputs: {input_cnt} outputs: {output_cnt}",
            f" catalysts: {catalyst_cnt} regulated: {regulatedby_cnt} inputs: {input_types} outputs: {output_types}",
            f" catalysts: {catalyst_types} regulated: {regulatedby_types}",
        )


def parse_complex_for_definedset(dbid, defined_sets: dict):

    doc = db.get_entity(dbid)
    # print("parse_complex", dbid, doc["schemaClass"])
    if "hasComponent" in doc:
        members = doc["hasComponent"]
    elif "hasMember" in doc:
        members = doc["hasMember"]
    else:
        members = []

    if doc["schemaClass"] == "DefinedSet":
        defined_sets[dbid] = len(members)

    for c in members:
        # print("C", c)
        if isinstance(c, int):
            c = db.get_entity(c)
        if c["schemaClass"] in ["DefinedSet", "Complex"]:
            parse_complex_for_definedset(c["dbId"], defined_sets)


def check_components(components: list):
    defined_sets = {}
    for c in components:
        if isinstance(c, int):
            c = db.get_entity(c)
        if c["schemaClass"] == "CatalystActivity":
            c = c["physicalEntity"]
            if isinstance(c, int):
                c = db.get_entity(c)
            if c["schemaClass"] in ["Complex", "DefinedSet"]:
                dbid = c["dbId"]
                parse_complex_for_definedset(dbid, defined_sets)

        elif c["schemaClass"] in ["Complex", "DefinedSet"]:
            dbid = c["dbId"]
            parse_complex_for_definedset(dbid, defined_sets)

    return defined_sets


def check_definedsets():
    reaction_ids = to_bel.collect_all_reactions()
    # reaction_ids = ["421007"]
    stats = {}
    count = 0
    for rid in reaction_ids:
        doc = db.get_entity(rid)
        dbid = doc["dbId"]
        schema_class = doc["schemaClass"]
        category = doc["category"]
        input_cnt = len(doc.get("input", []))
        output_cnt = len(doc.get("output", []))

        input_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("input", [])
            if isinstance(t, dict)
        ]
        output_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("output", [])
            if isinstance(t, dict)
        ]
        catalyst_cnt = len(doc.get("catalystActivity", []))
        catalyst_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("catalystActivity", [])
            if isinstance(t, dict)
        ]
        regulatedby_cnt = len(doc.get("regulatedBy", []))
        regulatedby_types = [
            f'{t.get("dbId")}:{t.get("schemaClass")}'
            for t in doc.get("regulatedBy", [])
            if isinstance(t, dict)
        ]

        input_definedsets = check_components(doc.get("input", []))
        output_definedsets = check_components(doc.get("output", []))
        catalyst_definedsets = check_components(doc.get("catalystActivity", []))

        if (
            len(input_definedsets.keys())
            + len(output_definedsets.keys())
            + len(catalyst_definedsets.keys())
            > 3
        ):
            print(f"\nRID: {dbid}   -- more than 3 definedsets in reaction")
            # else:
            #     print(f"\nRID: {dbid}")

            print(
                f"  class: {schema_class} category: {category} inputs: {input_cnt} outputs: {output_cnt}",
                f" catalysts: {catalyst_cnt} regulated: {regulatedby_cnt} inputs: {input_types} outputs: {output_types}",
                f" catalysts: {catalyst_types} regulated: {regulatedby_types}",
            )

            print(
                f"  DefinedSets Inputs: {json.dumps(input_definedsets)}  Outputs: {json.dumps(output_definedsets)}  Catalysts: {json.dumps(catalyst_definedsets)}"
            )
            if set(input_definedsets.keys()) != set(output_definedsets.keys()):
                print("MISMATCH - mismatched input/output definedsets")

        # count += 1
        # if count > 500:
        #     quit()


def review_complexes():
    reaction_ids = to_bel.collect_all_reactions()
    f = open("complexes.json", "w")
    stats = {}
    for rid in reaction_ids:
        print("Reaction ID", rid)
        doc = db.get_entity(rid)
        comps = to_bel.get_reaction_components(doc)
        keys = ["inputs", "outputs", "catalysts", "regulators"]
        for key in keys:
            if key in comps:
                for idx, item in enumerate(comps[key]):
                    if isinstance(item, to_bel.Complex):
                        flag = 0
                        for member in item.members:
                            if isinstance(member, to_bel.EntitySet):
                                flag += 1
                            elif isinstance(member, to_bel.Complex):
                                flag += 1

                        if flag > 2:
                            # Standard Library
                            import json

                            f.write(f"ReactionID: {rid}\n")
                            json.dump(item, f, cls=to_bel.CustomEncoder, indent=4)


def main():

    # TODO look for Traceback and Unable to process protein modification errors in logs

    # species_id = "9606"
    # pathways = collect.get_pathways(species_id=species_id)
    # reaction_ids = pathways_to_reactions(pathways=pathways)

    reaction_ids = []
    # reaction_ids = ["9615721"]

    to_bel.convert(reaction_ids, update_all=True, limit=1000000)  # Collect all human reactions
    sys.exit()

    # review_reactions()
    # quit()
    # review_complexes()
    # quit()

    # check_definedsets()
    # quit()

    # reaction_ids = ["9027627"]  # Catalyzed reaction
    # reaction_ids = ["445813"]  # Regulated by 445813
    # reaction_ids = ["2399988"]  # Disease annotation
    # # reaction_ids = ["2399988"]  # pmod() and var()

    reaction_ids = ["5655336"]  # has a defined/candidate 5655280 set as a subset
    # reaction_ids = ["8952044"]  # complexes of complexes with candidatesets
    # reaction_ids = ["8942101"]  # complexes of complexes of entity_sets

    to_bel.convert(reaction_ids)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted")
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
