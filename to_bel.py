import copy
import datetime
import glob
import itertools
import json
import re
from collections import Iterable
from typing import List

import db
import logging_setup
import settings
import structlog
from db import collect_all_reactions
from devtools import debug, pprint
from utils import clean_label, ts_convert

log = structlog.getLogger(__name__)

db_objects = db.get_db()
reactome_db = db_objects["reactome_db"]
reactome_db_name = db_objects["reactome_db_name"]
reactome_coll = db_objects["reactome_coll"]
reactome_coll_name = db_objects["reactome_coll_name"]
reactome_aql = db_objects["reactome_aql"]

sets = {}

unmatched_fh = open("unmatched.tsv", "w")
sets_fh = open("sets.json", "w")
errors_fh = open("errors.tsv", "w")

now = f"{datetime.datetime.utcnow().isoformat(timespec='milliseconds')}Z"

def quoting_entity_id(entity_id):
    """Quoting nsargs

    If needs quotes (only if it contains whitespace, comma or ')' ), make sure
        it is quoted, else don't add them.


    """
    quoted = re.findall(r'^"(.*)"$', entity_id)

    if re.search(
        r"[),\s]", entity_id
    ):  # quote only if it contains whitespace, comma or ')'
        if quoted:
            return entity_id
        else:
            return f'"{entity_id}"'
    else:
        if quoted:
            return quoted[0]
        else:
            return entity_id

class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (Complex, EntitySet, Regulator)):
            return obj.to_json()
        return json.JSONEncoder.default(self, obj)

class Entity:
    def __init__(
        self, namespace: str, id: str, label: str, stid: str, dbid: str, loc: str = ""
    ):
        self.namespace: str = namespace
        self.id: str = id
        self.label: str = label
        self.stid: str = stid
        self.dbid: str = dbid
        self.loc: str = loc

    def __str__(self):
        entity_id = quoting_entity_id(f"{self.id}!{self.label}")
        return f"{self.namespace}:{entity_id}"

    def __repr__(self):
        return self.__str__()


class Function:
    def __init__(
        self,
        function: str,
        parameters: List[Entity],
        stid: str,
        dbid: str,
        loc: str = "",
    ):
        self.function: str = function
        self.parameters: List[Entity] = parameters
        self.stid: str = stid
        self.dbid: str = dbid
        self.loc: str = loc

    def no_loc(self):
        parameters_list = []
        for parameter in self.parameters:
            if isinstance(parameter, str):
                parameters_list.append(parameter)
            else:
                parameters_list.append(str(parameter))
        return f'{self.function}({", ".join(parameters_list)})'

    def __str__(self):
        parameters_list = []
        for parameter in self.parameters:
            if isinstance(parameter, str):
                parameters_list.append(parameter)
            else:
                parameters_list.append(str(parameter))
        return f'{self.function}({", ".join(parameters_list)}, loc({self.loc}))'

    def __repr__(self):
        return self.__str__()


class Complex:
    def __init__(
        self, members: list, stid: str, dbid: str, label: str = "", loc: str = ""
    ):
        self.members: list = members
        self.stid: str = stid
        self.dbid: str = dbid
        self.label: str = label
        self.loc: str = loc

    def to_json(self):
        r = {
            "complex": {
                "members": self.members,
                "stid": self.stid,
                "dbid": self.dbid,
                "label": self.label,
                "loc": self.loc,
            }
        }
        return r

    def __str__(self):
        members_list = []
        for member in self.members:
            if isinstance(member, str):
                members_list.append(member)
            else:
                members_list.append(str(member))

        return f'complex({", ".join(members_list)}, loc({self.loc}))'

    def __repr__(self):
        return self.__str__()


class EntitySet:
    def __init__(
        self,
        members: list,
        stid: str,
        dbid: str,
        type: str,
        name: str = "",
        loc: str = "",
    ):
        self.members: list = members
        self.loc: str = loc
        self.stid: str = stid
        self.dbid: str = dbid
        self.name: str = name
        self.type: str = type

    def to_json(self):
        r = {
            "entity_set": {
                "members": self.members,
                "type": self.type,
                "name": self.name,
                "stid": self.stid,
                "dbid": self.dbid,
                "loc": self.loc,
            }
        }
        return r

    def __str__(self):
        members_list = []
        for member in self.members:
            if isinstance(member, str):
                members_list.append(member)
            else:
                members_list.append(str(member))

        return f'{self.type}({", ".join(members_list)}, loc({self.loc}))'

    def __repr__(self):
        return self.__str__()


class Regulator:
    def __init__(self, regulator, relation, dbid):
        self.regulator = regulator
        self.relation = relation
        self.dbid = dbid

    def to_json(self):
        r = {"regulator": self.regulator, "relation": self.relation, "dbid": self.dbid}
        return r

    def __str__(self):
        return f"Regulator: {self.regulator} Relation: {self.relation}"

    def __repr__(self):
        return self.__str__()


def get_creator(metadata, doc):

    if "authored" in doc:
        authored_id = doc["authored"][0]
        if isinstance(authored_id, dict):
            authored_id = authored_id["dbId"]

        authored_doc = db.get_entity(authored_id)

        if "author" in authored_doc:
            metadata["gd_creator"] = authored_doc["author"][0].get("displayName", "")
            metadata["gd_createTS"] = ts_convert(authored_doc["dateTime"])
            metadata["gd_creator_orcid"] = authored_doc["author"][0].get("orcidId", "")

    elif "created" in doc:
        matches = re.match(
            r"(.*?)\s+(\d{4,4}-\d{2,2}-\d{2,2}\s+\d{2,2}:\d{2,2}:\d{2,2})",
            doc["created"]["displayName"],
        )
        if matches:
            author = matches.group(1)
            timestamp = matches.group(2)
            timestamp = ts_convert(timestamp)
            metadata["gd_creator"] = author
            metadata["gd_createTS"] = timestamp
            return metadata

        matches = re.match(
            r"(.*?)\s+(\d{4,4}-\d{2,2}-\d{2,2})", doc["created"]["displayName"]
        )
        if matches:
            author = matches.group(1)
            timestamp = matches.group(2)
            timestamp = f"{timestamp} 00:00:00"
            timestamp = ts_convert(timestamp)
            metadata["gd_creator"] = author
            metadata["gd_createTS"] = timestamp
            return metadata

    return metadata


def get_species(species_id):
    doc = db.get_entity(species_id)
    return {"type": "Species", "id": f"TAX:{doc['taxId']}", "label": doc["displayName"]}


def get_disease(dbId):
    doc = db.get_entity(dbId)
    # TODO


def get_protein(doc) -> Entity:

    stid = doc["stId"]
    dbid = doc["dbId"]
    loc = ""

    try:
        namespace = ""
        if doc["referenceEntity"]["databaseName"] == "UniProt":
            namespace = "SP"
        else:
            log.info(
                f"Missing protein databaseName {doc['referenceEntity']['databaseName']}"
            )
            namespace = "TBD"

        label = clean_label(doc["referenceEntity"]["name"][0])
        id_: str = doc["referenceEntity"]["identifier"]

        loc = ""
        if doc.get("compartment", False):  # TODO - can there be multiple compartments?
            if doc["compartment"][0]["databaseName"] == "GO":
                loc_namespace = "GO"
            else:
                log.info(
                    f'Missing location databaseName {doc["compartment"][0]["databaseName"]}'
                )
                loc_namespace = "TBD"

            loc_id = doc["compartment"][0]["accession"]
            loc_label = clean_label(doc["compartment"][0]["displayName"])
            loc = f"{loc_namespace}:{loc_id}!{loc_label}"

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\tProtein\t{str(e)}\n")
        namespace = "REACTOME"
        id_ = dbid
        label = ""

    protein = Entity(namespace, id_, label, stid, dbid, loc=loc)
    p = Function("p", [protein], stid, dbid, loc=loc)

    # Debug
    debug(p)

    return p


def get_catalyst(doc) -> str:

    try:
        r = process_component(doc["physicalEntity"]["dbId"])

    except Exception as e:
        errors_fh.write(
            f"{doc['dbId']}\tCatalyst\tNot able to process physical entity for catalyst: {str(e)}\n"
        )
        r = ""

    return r


def get_gt(doc) -> str:
    """Get className Genes and Transcripts - seems like it's used for unknown proteins"""

    log.debug(f"GT className: {doc['dbId']}")

    prefix = "TBD"
    name = doc["name"][0]
    if "referenceEntity" in doc:
        name = clean_label(doc["referenceEntity"]["name"][0])
    else:
        name = clean_label(doc["name"][0])

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(
                f'Missing location databaseName {doc["compartment"][0]["databaseName"]}'
            )
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    if loc:
        protein = f"p({prefix}:{name}, loc({loc}))"
    else:
        protein = f"p({prefix}:{name})"

    return protein


def get_compound(doc: dict) -> str:

    if doc["referenceEntity"]["databaseName"] == "ChEBI":
        prefix = "CHEBI"
    elif doc["referenceEntity"]["databaseName"] == "IUPHAR":
        prefix = "IUPHAR"
    else:
        log.info(
            f"Missing compound databaseName {doc['referenceEntity']['databaseName']}"
        )
        prefix = "TBD"

    name = clean_label(doc["referenceEntity"]["name"][0])
    _id = doc["referenceEntity"]["identifier"]

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(
                f'Missing location databaseName {doc["compartment"][0]["databaseName"]}'
            )
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    if loc:
        compound = f"a({prefix}:{_id}!{name}, loc({loc}))"
    else:
        compound = f"a({prefix}:{_id}!{name})"

    return compound


def get_complex(doc: dict) -> str:
    """Get complex"""

    components = set()
    for component in doc["hasComponent"]:
        if not isinstance(component, dict):
            continue

        r = process_component(component["dbId"])
        components.add(r)

    stid = doc["stId"]
    dbid = doc["dbId"]

    name = doc.get("name", [""])[0]
    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(
                f'Missing location databaseName {doc["compartment"][0]["databaseName"]}'
            )
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    c = Complex(list(components), name=name, dbid=dbid, stid=stid, loc=loc)

    return c


def get_set(doc: dict) -> dict:

    json.dump(doc, sets_fh, indent=4)
    # log.info("Set", doc=doc)

    name = doc["name"][0]
    name = clean_label(name)
    protein = f"REACTOME:{name}"
    members = {}

    if "hasMember" in doc:
        for member in doc["hasMember"]:
            r = process_component(member["dbId"])
            members[r] = 1

    if "hasCandidate" in doc:
        for member in doc["hasCandidate"]:
            r = process_component(member["dbId"])
            members[r] = 1

    stid = doc["stId"]
    dbid = doc["dbId"]
    name = doc.get("name", [""])[0]

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(
                f'Missing location databaseName {doc["compartment"][0]["databaseName"]}'
            )
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    s = EntitySet(
        list(members.keys()),
        stid=stid,
        dbid=dbid,
        name=name,
        type=doc["schemaClass"],
        loc=loc,
    )

    return s


def get_regulator(doc) -> str:
    """Get regulation"""

    relation = "regulates"
    try:
        if isinstance(doc["regulator"], int):
            rdbid = doc["regulator"]
        else:
            rdbid = doc["regulator"]["dbId"]
        regulator = process_component(rdbid)

        if doc["className"] in [
            "PositiveGeneExpressionRegulation",
            "PositiveRegulation",
        ]:
            relation = "increases"
        elif doc["className"] in [
            "NegativeGeneExpressionRegulation",
            "NegativeRegulation",
        ]:
            relation = "decreases"
        else:
            log.info(f"Unknown regulation relationship: {doc['dbId']}")

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\tRegulator\t{str(e)}\n")
        regulator = f'Cannot process Regulator for {doc["regulator"]["dbId"]}'

    # print("regulator", regulator, "relation", relation)

    dbid = doc["dbId"]
    r = Regulator(regulator, relation, dbid)
    return r


def get_gene(doc) -> str:
    """Get gene"""

    try:
        prefix = doc["referenceEntity"]["databaseName"]
        _id = doc["referenceEntity"]["identifier"]
        name = doc["referenceEntity"]["geneName"][0]
        name = clean_label(name)

        gene = f"g({prefix}:{_id}!{name})"

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\tGene\t{str(e)}\n")
        gene = f"g(REACTOME:{doc['dbId']})"

    return gene


def get_rna(doc) -> str:
    """Get RNA"""

    try:
        prefix = doc["referenceEntity"]["databaseName"]
        _id = doc["referenceEntity"]["identifier"]
        name = doc["referenceEntity"]["geneName"][0]
        name = clean_label(name)

        rna = f"r({prefix}:{_id}!{name})"

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\RNA\t{str(e)}\n")
        rna = f"r(REACTOME:{doc['dbId']})"

    return rna


def get_polymer(doc) -> str:
    """Get Polymer"""

    prefix = "REACTOME"
    stid = doc["stId"]
    name = doc["name"][0]
    name = clean_label(name)

    polymer = f"r({prefix}:{stid}!{name})"

    return polymer


def get_other_entity(doc) -> str:
    """Get OtherEntity"""

    prefix = "REACTOME"
    name = doc["name"][0]
    stid = doc["stId"]
    name = clean_label(name)

    other_entity = f"FNTBD({prefix}:{stid}!{name})"

    return other_entity


def get_requirement(doc):
    if isinstance(doc["regulator"], int):
        dbid = doc["regulator"]
    else:
        dbid = doc["regulator"]["dbId"]
    doc = db.get_entity(dbid)

    if doc["className"] == "Complex":
        return get_complex(doc)
    elif doc["className"] == "Chemical Compound":
        return get_compound(doc)


def process_component(dbid: str) -> dict:
    """Process reaction component: input, output, catalyst or regulator"""

    doc = db.get_entity(dbid)

    if doc["className"] == "Protein":
        return get_protein(doc)
    elif doc["className"] == "Complex":
        return get_complex(doc)
    elif doc["className"] == "Requirement":
        return get_requirement(doc)
    elif doc["className"] in ["Chemical Compound", "ChemicalDrug", "ProteinDrug"]:
        return get_compound(doc)
    elif doc["className"] == "Set":
        return get_set(doc)
    elif doc["className"] == "DNA Sequence":
        return get_gene(doc)
    elif doc["className"] == "RNA Sequence":
        return get_rna(doc)
    elif doc["className"] == "Genes and Transcripts":
        return get_gt(doc)
    elif doc["className"] == "Polymer":
        return get_polymer(doc)
    elif doc["className"] == "OtherEntity":
        return get_other_entity(doc)
    elif doc["className"] == "CatalystActivity":
        return get_catalyst(doc)
    elif doc["className"] in [
        "PositiveGeneExpressionRegulation",
        "NegativeGeneExpressionRegulation",
        "NegativeRegulation",
        "PositiveRegulation",
    ]:
        return get_regulator(doc)
    else:
        log.info(f'Unmatched class name {doc["className"]} for dbId: {dbid}')
        unmatched_fh.write(f"{dbid}\t{doc['className']}\t{doc['displayName']}\n")
        return ""


# reaction types: [
#   "transition",
#   "omitted",
#   "binding",
#   "dissociation",
#   "uncertain"
# ]
def get_reaction_components(doc: dict) -> dict:

    doc = db.get_entity(doc["dbId"])
    desc = doc["displayName"]
    reaction_type = doc["category"]

    inputs = []
    outputs = []
    catalysts = []
    regulators = []

    # inputs
    for input in doc.get("input", []):
        if isinstance(input, int):
            dbid = input
        else:
            dbid = input["dbId"]
        r = process_component(dbid)
        inputs.append(r)

    for output in doc.get("output", []):
        if isinstance(output, int):
            dbid = output
        else:
            dbid = output["dbId"]
        r = process_component(dbid)
        outputs.append(r)

    for regulator in doc.get("regulatedBy", []):
        if isinstance(regulator, int):
            dbid = regulator
        else:
            dbid = regulator["dbId"]
        r = process_component(dbid)
        regulators.append(r)

    for catalyst in doc.get("catalystActivity", []):
        if isinstance(catalyst, int):
            dbid = catalyst
        else:
            dbid = catalyst["dbId"]
        r = process_component(dbid)
        catalysts.append(r)

    inputs = list(set(inputs))
    outputs = list(set(outputs))
    catalysts = list(set(catalysts))
    regulators = list(set(regulators))

    return {
        "dbId": doc["dbId"],
        "description": desc,
        "reaction_type": reaction_type,
        "inputs": inputs,
        "outputs": outputs,
        "catalysts": catalysts,
        "regulators": regulators,
    }


# def flatten_obj(obj, set_members: list, complex_members: list):
#     """Flatten complexes and entity to non-redundant set"""

#     for member in obj.members:
#         if isinstance(member, (Complex, EntitySet)):
#             sm, cm = flatten_obj(member, set_members, complex_members)
#             set_members.extend(sm)
#             complex_members.extend(cm)

#     return set_members, complex_members


# def flatten(items):
#     """Yield items from any nested iterable; see Reference."""
#     for x in items:
#         if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
#             for sub_x in flatten(x):
#                 yield sub_x
#         else:
#             yield x


# def sort_complex(members):
#     """Sort complex members and remove redundancies"""

#     members = sorted(list(flatten(members)))

#     return members


# # def flatten_set(entity_set):
# #     """Flatten entity_set object"""

# #     global sets

# #     members = []
# #     for member in entity_set.members:
# #         if isinstance(member, EntitySet):
# #             members.append(flatten_set(member))
# #         elif isinstance(member, Complex):
# #             members.append(flatten_complex(member))
# #         else:
# #             members.append(member)

# #     return members


def get_function_type(dbid, members):
    """Get entity set function type"""
    function_types = []
    for member in members:
        matches = re.match("^(\w+)\(", member)
        prefix = matches.group(1)
        function_types.append(prefix)

    function_types = list(set(function_types))
    if len(function_types) == 1:
        return function_types[0]
    elif len(function_types) > 0 and "complex" in function_types:
        log.warning(f"Multiple function types for Entity Set {dbid}: {function_types}")
        return "complex"
    else:
        log.warning(f"Multiple function types for Entity Set {dbid}: {function_types}")

    return "MISSING"


# def format_objects(comps):
#     """Format complex and entity_sets into """

#     global sets

#     keys = ["inputs", "outputs", "catalysts", "regulators"]
#     for key in keys:
#         print("Key", key)
#         if key in comps:
#             for idx, item in enumerate(comps[key]):
#                 print("Item", item)
#                 if isinstance(item, EntitySet):
#                     log.info("Formatting entity sets")
#                     name = item.name
#                     stid = item.stid
#                     dbid = item.dbid
#                     set_members, complex_members = flatten_obj(item, [], [])
#                     members = list(set(flatten(set_members)))
#                     comps[key][idx] = set_members
#                     # function_type = get_function_type(dbid, members)
#                     # bel_fn = f"{function_type}(reactome:{stid})"
#                     # comps[key][idx] = bel_fn
#                     # sets[bel_fn] = []
#                     # for member in members:
#                     #     sets[bel_fn].append({"subject": bel_fn, "relation": "hasMember", "object": member})

#                 elif isinstance(item, Complex):
#                     log.info("Formatting complex")
#                     # members = sort_complex(flatten_obj(item))
#                     set_members, complex_members = flatten_obj(item, [], [])
#                     members = ", ".join(complex_members)
#                     comps[key][idx] = f"complex({members})"

#                 elif isinstance(item, Regulator):
#                     if isinstance(item.regulator, Complex):
#                         log.info("Formatting complex")
#                         # members = sort_complex(flatten_obj(item.regulator))
#                         set_members, complex_members = flatten_obj(item, [], [])
#                         members = ", ".join(complex_members)
#                         item.regulator = f"complex({members})"
#                         comps[key][idx] = item

#     return comps


def process_candidate_set(candidateset):

    # print("Processing CandidateSet")
    dbid = candidateset.dbid
    stid = candidateset.stid
    name = clean_label(candidateset.name)
    fn = get_function_type(dbid, candidateset.members)
    sets[dbid] = []
    for member in candidateset.members:
        sets[dbid].append(
            {
                "subject": f"{fn}(REACTOME:{stid}!{name})",
                "relation": "hasMember",
                "object": member,
            }
        )

    # debug(sets[dbid])
    return f"{fn}(REACTOME:{stid}!{name})"


def format_complex(complex, complex_list: list, level: int = 0):
    """Format complex"""

    global sets

    # print("Start Complex", complex)

    normal_members = []
    set_members = []

    for member in complex.members:
        if isinstance(member, Complex):
            # print(f"{' ' * level* 2}Complex member {member.dbid}")
            debug(member.members)
            level += 1
            complex_list.extend(format_complex(member, complex_list, level))
        elif isinstance(member, EntitySet):
            if member.type == "CandidateSet":
                candidateset = process_candidate_set(member)
                # print(f"{' ' * level * 2}CandidateSet {member.dbid}", candidateset)
                debug(candidateset)
                normal_members.append(candidateset)

            else:
                # print(f"{' ' * level * 2}DefinedSet {member.dbid}")
                debug(member.members)
                set_members.append(member.members)
        else:
            # print(f"{' ' * level * 2}Member")
            debug(member)
            normal_members.append(member)

    # debug(normal_members)
    # debug(set_members)

    temp_list = []
    enumerated_set_members = list(itertools.product(*set_members))
    if enumerated_set_members:
        for esm in enumerated_set_members:
            temp_list.append(list(esm) + normal_members)
    else:
        temp_list = normal_members

    # debug(temp_list)
    # debug(complex_list)

    if complex_list and temp_list:
        new_complex_list = []
        for cl in complex_list:
            for tl in temp_list:
                new_complex_list.append(cl + tl)
        return new_complex_list
    elif complex_list:
        return complex_list
    elif temp_list:
        return temp_list


def format_comps(comps):

    global sets
    keys = ["catalysts", "regulators", "inputs", "outputs"]
    for key in keys:
        print(f"\n\n#############  Key {key}  ################")
        if key in comps:
            for idx, item in enumerate(comps[key]):
                complexes = []
                print("Item", item)
                if isinstance(item, Complex):
                    complex_members = format_complex(item, [], 0)
                    for cm in complex_members:
                        complexes.append(f"complex({', '.join(sorted(list(set(cm))))})")
                    comps[key][idx] = list(set(complexes))
                    # debug(complexes)
                if isinstance(item, EntitySet):
                    if item.type == "CandidateSet":
                        dbid = item.dbid
                        stid = item.stid
                        name = clean_label(item.name)

                        comps[key][idx] = f"{fn}(REACTOME:{stid}!{name})"
                    else:
                        comps[key][idx] = item.members
                if isinstance(item, Regulator):
                    if isinstance(item.regulator, Complex):
                        complex_members = format_complex(item.regulator, [], 0)
                        for cm in complex_members:
                            complexes.append(
                                f"complex({', '.join(sorted(list(set(cm))))})"
                            )
                        comps[key][idx].regulator = list(set(complexes))

    # debug(comps)


def create_assertions(doc):

    comps = get_reaction_components(doc)
    debug(comps)
    comps = format_comps(comps)
    debug(comps)
    assertions = []

    return assertions


def process_reaction(reaction_id):

    processing_reactions_flag = False

    doc = db.get_entity(reaction_id)
    if doc["schemaClass"] == "FailedReaction":
        # TODO - handle Failed Reaction - e.g. 3229118 (http://localhost:18529/_db/reactome_db/_admin/aardvark/index.html#collection/reactome/3229118)
        return {}

    src_url = f"https://reactome.org/content/detail/{doc['stId']}"
    metadata = {
        "source_url": src_url,
        "source": "Reactome",
        "license": "CC0",
        "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
        "gd_updateTS": now,
    }

    log.info("Reaction", dbid={doc["dbId"]})

    metadata = get_creator(metadata, doc)

    nanopub = {"assertions": [], "annotations": []}
    annotations = []
    nanopub["id"] = f"Reactome_{doc['stId']}"

    nanopub["evidence"] = doc["summation"][0]["text"]
    nanopub["citation"] = {"uri": src_url}

    if "species" in doc:
        species_id = doc["species"][0]
        if isinstance(species_id, dict):
            species_id = species_id["dbId"]

        annotations.append(
            get_species(species_id)
        )  # TODO - can there be multiple species?

    if "disease" in doc:
        annotations.append(
            get_disease(doc["disease"][0]["dbId"])
        )  # TODO - can there be multiple diseases?

    nanopub["assertions"] = create_assertions(doc)
    nanopub["annotations"] = copy.deepcopy(annotations)
    nanopub["metadata"] = copy.deepcopy(metadata)

    return nanopub


def convert(reaction_ids: list = []):
    """Convert reactome reactions to BEL

    If reaction_ids is empty - it will process all Reactions in database
    """

    if not reaction_ids:
        reaction_ids = collect_all_reactions()

    f = open("reactome_nanopubs.jsonl", "w")
    try:
        for reaction_id in reaction_ids:
            nanopub = process_reaction(reaction_id)

            print("Nanopub:\n", json.dumps(nanopub, indent=4))

            if not nanopub:
                continue

            # reactions.append(nanopub)
            f.write(f"{json.dumps(nanopub, cls=CustomEncoder)}\n")

    except Exception as e:
        log.exception("Converting reactions", error=e)
    finally:
        f.close()
