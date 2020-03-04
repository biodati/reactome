import copy
import datetime
import glob
import itertools
import json
import re
from collections.abc import Iterable
from typing import List, Optional

import xxhash
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

entity_sets = {}  # named families to write into BEL statemens (isA relations)
# complexes = {}  # named complexes to write into BEL statements (hasComponent relations)

now = f"{datetime.datetime.utcnow().isoformat(timespec='milliseconds')}Z"


def quote_entity(entity_str):
    """Quoting nsargs

    If needs quotes (only if it contains whitespace, comma or ')' ), make sure
        it is quoted, else don't add them.


    """
    entity_str = str(entity_str)

    quoted = re.findall(r'^"(.*)"$', entity_str)

    if re.search(r"[),\s]", entity_str):  # quote only if it contains whitespace, comma or ')'
        if quoted:
            return entity_str
        else:
            return f'"{entity_str}"'
    else:
        if quoted:
            return quoted[0]
        else:
            return entity_str


class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (Function, Entity, Complex, EntitySet, Regulator, Catalyst, PMod)):
            return obj.to_json()
        return json.JSONEncoder.default(self, obj)


class Entity:
    def __init__(
        self,
        namespace: str,
        id: str,
        stid: str,
        stid_version: str,
        dbid: str,
        label: str = "",
        loc: str = "",
    ):
        self.namespace: str = namespace
        self.id: str = id
        self.label: str = label
        self.stid: str = stid
        self.stid_version: str = stid_version
        self.dbid: str = dbid
        self.loc: str = loc
        self.type = "Entity"

    def to_json(self):
        r = {
            "entity": {
                "namespace": self.namespace,
                "id": self.id,
                "stid_version": self.stid_version,
                "dbid": self.dbid,
                "label": self.label,
                "loc": self.loc,
            }
        }
        return r

    def to_bel(self, strip_locations=False):

        if self.id == self.label or not self.label:
            entity_id = quote_entity(self.id)
            return f"{self.namespace}:{entity_id}"
        else:
            entity_id = quote_entity(self.id)
            entity_label = quote_entity(self.label)
            return f"{self.namespace}:{entity_id}!{entity_label}"

    def __str__(self):
        return self.to_bel()

    def __repr__(self):
        return self.__str__()


class Function:
    def __init__(
        self,
        function: str,
        parameters: List[Entity],
        stid: str,
        stid_version: str,
        dbid: str,
        mods: list = [],
        loc: str = "",
    ):
        self.function: str = function
        self.parameters: List[Entity] = parameters
        self.stid: str = stid
        self.stid_version: str = stid_version
        self.dbid: str = dbid
        self.loc: str = loc
        self.type = "Function"
        self.mods = mods

    def to_json(self):
        r = {
            "Function": {
                "function": self.function,
                "parameters": self.parameters,
                "stid": self.stid,
                "stid_version": self.stid_version,
                "dbid": self.dbid,
                "loc": self.loc,
                "mods": self.mods,
            }
        }
        return r

    def to_bel(self, strip_locations=False):
        parameters_list = []
        for parameter in self.parameters:
            if isinstance(parameter, str):
                parameters_list.append(parameter)
            else:
                parameters_list.append(parameter.to_bel(strip_locations=strip_locations))        

        location = ""
        if not strip_locations and self.loc:
            location = f", loc({self.loc})"

        modifications = []
        if self.mods:
            for mod in self.mods:
                modifications.append(mod.to_bel())

            modifications = f', {", ".join(modifications)}'
        else:
            modifications = ""
            
        bel = f'{self.function}({", ".join(parameters_list)}{modifications}{location})'

        log.debug("Function to_bel", bel=bel, strip_locations=strip_locations)

        return bel

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

    def __lt__(self, other):
        return self.__repr__() < other.__repr__()

    def __le__(self, other):
        return self.__repr__() <= other.__repr__()

    def __gt__(self, other):
        return self.__repr__() > other.__repr__()

    def __ge__(self, other):
        return self.__repr__() >= other.__repr__()

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def __ne__(self, other):
        return self.__repr__() != other.__repr__()

    def  __hash__(self):
        return int(xxhash.xxh32(str(self.__repr__())).intdigest())


class PMod:
    def __init__(self, modification, dbid, residue: Optional[str] = None, coordinate: Optional[int] = None):
        self.modification = modification
        self.residue = residue
        self.coordinate = coordinate
        self.dbid = dbid
    
    def to_bel(self, strip_locations=False):

        if self.modification and self.residue and self.coordinate:
            return f"pmod({self.modification}, {self.residue}, {self.coordinate})"
        elif self.modification and self.residue:
            return f"pmod({self.modification}, {self.residue})"
        else:
            return f"pmod({self.modification})"
    
    def mod_only_to_bel(self):
        """Only provide the modification - used for EntitySets"""

        return f"pmod({self.modification})"

    def to_json(self):
        r = {
            "pmod": {
                "modification": self.modification,
                "residue": self.residue,
                "coordinate": self.coordinate,
                "dbid": self.dbid,
            }
        }
        return r

    def __str__(self):
        return self.to_bel()

    def __repr__(self):
        return self.__str__()

class Complex:
    def __init__(
        self, members: list, stid: str, stid_version: str, dbid: str, label: str = "", loc: str = ""
    ):
        self.members: list = members
        self.function = "complex"
        self.stid: str = stid
        self.stid_version: str = stid_version
        self.dbid: str = dbid
        self.label: str = label
        self.type = "Complex"
        self.loc: str = loc

        self.get_location()

    def to_json(self):
        r = {
            "complex": {
                "members": self.members,
                "stid": self.stid,
                "stid_version": self.stid_version,
                "dbid": self.dbid,
                "label": self.label,
                "loc": self.loc,
            }
        }
        return r

    def get_members(self):
        """Recursively return Function members of complexes"""
        pass

    def get_location(self):

        locations = set()
        for member in self.members:
            if member.type == "Function":
                locations.add(member.loc)

        log.debug("Locations", locations)

        if len(locations) == 1:
            self.loc = list(locations)[0]
        else:
            log.warning("Complex has mis-matched locations", complex=self.__repr__)
            self.loc = ""
    
    def get_members(self):
        """Get Complex members to flatten out the complex members"""

        members = set()
        for member in self.members:
            if isinstance(member, Complex):
                for c in member.get_members():
                    members.add(c)
            else:
                members.add(member)
        
        return list(members)
        
    def to_bel(self, strip_locations=True):
        strip_locations = True  # override this variable

        self.get_location()

        members_strings = []
        for member in sorted(self.get_members()):
            members_strings.append(member.to_bel(strip_locations=strip_locations))

        return f"complex({', '.join(members_strings)}, loc({self.loc})"

    def __str__(self):
        return self.to_bel()

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return self.__repr__() < other.__repr__()
    
    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def  __hash__(self):
        return int(xxhash.xxh32(str(self.__repr__())).intdigest())
        

class EntitySet:
    def __init__(
        self,
        id: str,
        members: list,
        stid: str,
        stid_version: str,
        dbid: str,
        type_: str,
        namespace: str = "REACTOME",
        label: str = "",
        loc: str = "",
        function: str = "",
    ):
        self.namespace = namespace
        self.id = id
        self.members: list = members
        self.stid: str = stid
        self.stid_version: str = stid_version
        self.dbid: str = dbid
        self.label: str = label
        self.type: str = type_

        self.get_location()
        self.get_function()
        self.get_modifications()

    def to_json(self):
        r = {
            "entity_set": {
                "namespace": self.namespace,
                "id": self.id,
                "members": self.members,
                "type": self.type,
                "label": self.label,
                "stid": self.stid,
                "stid_version": self.stid_version,
                "dbid": self.dbid,
                "loc": self.loc,
            }
        }
        return r

    def get_function(self):
        functions = set()
        for member in self.members:
            if isinstance(member, Function):
                functions.add(member.function)
            else:
                log.warning(f"EntitySet {self.dbid} has non-function member: {str(member)}")
        if len(functions) == 1:
            self.function = list(functions)[0]
        else:
            self.function = "MISSING"
            log.warning(f"EntitySet {self.dbid} has multiple or no functions {str(list(functions))}")

    def get_location(self):
        locations = set()
        for member in self.members:
            if isinstance(member, Function):
                locations.add(member.loc)
            else:
                log.warning(f"EntitySet {self.dbid} has no location: {str(member)}")

        if len(locations) == 1:
            self.loc = list(locations)[0]
        else:
            log.warning(f"EntitySet {self.dbid} has multiple or no locations {str(list(locations))}")
            self.loc = ""

    def get_modifications(self):
        mods = set()
        for member in self.members:
            if isinstance(member, Function):
                for mod in member.mods:
                    mods.add(mod)
            else:
                log.warning(f"EntitySet {self.dbid} has non-function member: {str(member)}")

        if len(mods) > 0:
            # TODO - simplify mods
            self.mods = list(mods)
        else:
            # log.warning(f"EntitySet {self.dbid} has multiple or no locations {str(list(mods))}")
            self.mods = []

    def to_bel(self, strip_locations=False):

        location = ""
        if not strip_locations and self.loc:
            location = f", loc({self.loc})"

        modifications = ""
        if self.mods:
            modifications = list(set([mod.mod_only_to_bel() for mod in self.mods]))
            modifications = ", ".join(modifications)
            modifications = f", {modifications}"

        return f"{self.function}(REACTOME:{self.stid_version}!{quote_entity(self.label)}{modifications}{location})"

    def list_members(self):
        members_list = []
        for member in self.members:
            if isinstance(member, str):
                members_list.append(member)
            else:
                members_list.append(str(member))

        return f'{self.type}({", ".join(members_list)}, loc({self.loc}))'


    def __str__(self):
        return f"REACTOME:{self.stid_version}!{quote_entity(self.label)} - {self.type}"


    def __repr__(self):
        return self.__str__()


class Regulator:
    def __init__(
        self, regulator, relation, dbid, stid, stid_version: str,
    ):
        self.regulator = regulator
        self.relation = relation
        self.dbid = dbid
        self.stid = stid
        self.stid_version: str = stid_version
        self.type = "Regulator"

    def to_json(self):
        r = {
            "regulator": self.regulator,
            "relation": self.relation,
            "dbid": self.dbid,
            "stid": self.stid,
            "stid_version": self.stid_version,
        }
        return r

    def to_bel(self, strip_locations=False):
        return self.__repr__()

    def __str__(self):
        return f"Regulator: {self.regulator} Relation: {self.relation}"

    def __repr__(self):
        return self.__str__()


class Catalyst:
    def __init__(
        self, catalyst, relation, dbid, stid, stid_version: str,
    ):
        self.catalyst = catalyst
        self.relation = relation
        self.dbid = dbid
        self.stid = stid
        self.stid_version: str = stid_version
        self.type = "Catalyst"

    def to_json(self):
        r = {
            "catalyst": self.catalyst,
            "relation": self.relation,
            "dbid": self.dbid,
            "stid": self.stid,
            "stid_version": self.stid_version,
        }
        return r

    def to_bel(self, strip_locations=False):
        return self.__repr__()

    def __str__(self):
        return f"Catalyst: {self.catalyst} Relation: {self.relation}"

    def __repr__(self):
        return self.__str__()


def get_creator(metadata, doc):

    if "authored" in doc:
        authored_id = doc["authored"][0]
        if isinstance(authored_id, dict):
            authored_id = authored_id["dbId"]

        authored_doc = db.get_entity(authored_id)

        if "author" in authored_doc:
            metadata["creator"] = authored_doc["author"][0].get("displayName", "")
            metadata["gd_createTS"] = ts_convert(authored_doc["dateTime"])
            metadata["creator_orcid"] = authored_doc["author"][0].get("orcidId", "")

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

        matches = re.match(r"(.*?)\s+(\d{4,4}-\d{2,2}-\d{2,2})", doc["created"]["displayName"])
        if matches:
            author = matches.group(1)
            timestamp = matches.group(2)
            timestamp = f"{timestamp} 00:00:00"
            timestamp = ts_convert(timestamp)
            metadata["creator"] = author
            metadata["gd_createTS"] = timestamp
            return metadata

    return metadata


def get_species(species_id):
    doc = db.get_entity(species_id)
    return {"type": "Species", "id": f"TAX:{doc['taxId']}", "label": doc["displayName"]}


def get_disease(dbId):
    doc = db.get_entity(dbId)
    # TODO


def process_mod(dbid):
    doc = db.get_entity(dbid)

    if "psiMod" in doc:
        if isinstance(doc["psiMod"], dict):
            psimods = [doc["psiMod"]]
        else:
            psimods = doc["psiMod"]

        for psimod in psimods:
            import json
            print('DumpVar:\n', json.dumps(psimod, indent=4))
            modification = ""
            residue = ""
            coordinate = psimod.get("coordinate", "")

            if psimod["identifier"] == "00046":
                modification = "Ph"
                residue = "Ser"

            elif psimod["identifier"] == "00047":
                modification = "Ph"
                residue = "Thr"

            else:
                modification = f'PSIMOD:{psimod["identifier"]}'
                
            return PMod(modification=modification, residue=residue, coordinate=coordinate, dbid=dbid)

        else:
            log.error(f"Unable to process modification: {dbid}")
            return PMod("Missing")


def get_protein(doc) -> Entity:

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    loc = ""
    mods = []

    try:
        namespace = ""
        if doc["referenceEntity"]["databaseName"] == "UniProt":
            namespace = "SP"
        else:
            log.info(f"Missing protein databaseName {doc['referenceEntity']['databaseName']}")
            namespace = "TBD"

        label = clean_label(doc["referenceEntity"]["name"][0])
        id_: str = doc["referenceEntity"]["identifier"]

        loc = ""
        if doc.get("compartment", False):  # TODO - can there be multiple compartments?
            if doc["compartment"][0]["databaseName"] == "GO":
                loc_namespace = "GO"
            else:
                log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
                loc_namespace = "TBD"

            loc_id = doc["compartment"][0]["accession"]
            loc_label = clean_label(doc["compartment"][0]["displayName"])
            loc = f"{loc_namespace}:{loc_id}!{loc_label}"

        if "hasModifiedResidue" in doc:
            for modification in doc["hasModifiedResidue"]:
                mods.append(process_mod(modification["dbId"]))

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\tProtein\t{str(e)}\n")
        log.exception(f"{doc['dbId']}\tProtein\tError: {str(e)}\n")
        namespace = "REACTOME"
        id_ = dbid
        label = ""

    protein = Entity(namespace=namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid, loc=loc)
    p = Function(function="p", parameters=[protein], mods=mods, stid=stid, stid_version=stid_version, dbid=dbid, loc=loc)

    return p


def get_catalyst(doc) -> str:
    """Get catalyst"""

    # print("Dump catalyst doc:\n", json.dumps(doc, indent=4))

    try:
        dbid = doc["physicalEntity"]["dbId"]
        stid = doc["physicalEntity"]["stId"]
        stid_version = doc["physicalEntity"]["stIdVersion"]

        catalyst = process_component(dbid)
        if catalyst.function != "a":
            catalyst = Function(function="activity", parameters=[catalyst], stid=stid, stid_version=stid_version, dbid=dbid)
        

        # relation
        relation = "directlyIncreases"

        r = Catalyst(catalyst, relation=relation, dbid=dbid, stid=stid, stid_version=stid_version)

    except Exception as e:
        errors_fh.write(f"{doc['physicalEntity']['dbId']}\tCatalyst\t{str(e)}\n")
        log.exception(f"{doc['physicalEntity']['dbId']}\tCatalyst\tError: {str(e)}\n")
        catalyst = f"Cannot process Catalyst for {doc['physicalEntity']['dbId']}"
        r = Catalyst(catalyst, relation="", dbid="", stid="", stid_version="")

    return r


def get_gt(doc) -> str:
    """Get className Genes and Transcripts - seems like it's used for unknown proteins"""

    log.debug(f"GT className: {doc['dbId']}")

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    loc = ""

    namespace = "TBD"
    label = doc["name"][0]
    if "referenceEntity" in doc:
        label = clean_label(doc["referenceEntity"]["name"][0])
    else:
        label = clean_label(doc["name"][0])

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    protein = Entity(namespace, id=label, label=label, stid=stid, stid_version=stid_version, dbid=dbid, loc=loc)
    p = Function("p", [protein], stid, stid_version, dbid, loc=loc)

    return p


def get_compound(doc: dict) -> str:

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    loc = ""

    if doc["referenceEntity"]["databaseName"] == "ChEBI":
        namespace = "CHEBI"
    elif doc["referenceEntity"]["databaseName"] == "IUPHAR":
        namespace = "IUPHAR"
    else:
        log.info(f"Missing compound databaseName {doc['referenceEntity']['databaseName']}")
        namespace = "TBD"

    label = clean_label(doc["referenceEntity"]["name"][0])
    id_ = doc["referenceEntity"]["identifier"]

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    compound = Entity(namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid, loc=loc)
    c = Function("a", [compound], stid, stid_version, dbid, loc=loc)

    return c


def get_complex(doc: dict) -> str:
    """Get complex"""

    components = set()
    for component in doc["hasComponent"]:
        if isinstance(component, int):
            dbid = component
        else:
            dbid = component["dbId"]

        r = process_component(dbid)

        components.add(r)

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]

    label = doc.get("name", [""])[0]
    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    c = Complex(
        list(components), label=label, dbid=dbid, stid=stid, stid_version=stid_version, loc=loc
    )

    return c


def get_set(doc: dict) -> dict:

    json.dump(doc, sets_fh, indent=4)
    # log.info("Set", doc=doc)

    label = doc["name"][0]
    label = clean_label(label)
    members = set()

    if "hasMember" in doc:
        for member in doc["hasMember"]:
            r = process_component(member["dbId"])
            members.add(r)

    if "hasCandidate" in doc:
        for member in doc["hasCandidate"]:
            r = process_component(member["dbId"])
            members.add(r)

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    label = doc.get("name", [""])[0]

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        loc = f"{loc_prefix}:{loc_id}!{loc_name}"

    s = EntitySet(
        members=list(members),
        namespace="REACTOME",
        id=stid_version,
        stid=stid,
        stid_version=stid_version,
        dbid=dbid,
        label=label,
        type_=doc["schemaClass"],
        loc=loc,
    )

    # Adding 
    entity_sets[dbid] = s

    return s


def get_regulator(doc) -> str:
    """Get regulation"""

    # print("DumpVar:\n", json.dumps(doc, indent=4))

    if 'regulator' not in doc:
        log.error("Bad regulator doc - missing regulator attribute", dbid=doc["dbId"])
        return Regulator("", relation="", dbid="", stid="", stid_version="")

    try:
        dbid = doc["regulator"]["dbId"]
        stid = doc["regulator"]["stId"]
        stid_version = doc["regulator"]["stIdVersion"]

        regulator = process_component(dbid)

        # Wrap as activity
        if regulator.function != "a":
            regulator = Function(function="activity", parameters=[regulator], stid=stid, stid_version=stid_version, dbid=dbid)
        
        # relation
        if doc["className"] in [
            "PositiveGeneExpressionRegulation",
            "PositiveRegulation",
        ]:
            relation = "directlyIncreases"

        elif doc["className"] in [
            "NegativeGeneExpressionRegulation",
            "NegativeRegulation",
        ]:
            relation = "directlyDecreases"
        else:
            log.info(f"Unknown regulation relationship: {doc['dbId']}")
            relation = "regulates"

        r = Regulator(regulator, relation=relation, dbid=dbid, stid=stid, stid_version=stid_version)

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\tRegulator\t{str(e)}\n")
        log.exception(f"{doc['dbId']}\tRegulator\tError: {str(e)}\n")
        regulator = f'Cannot process Regulator for {doc["regulator"]["dbId"]}'
        r = Regulator(regulator, relation="", dbid="", stid="", stid_version="")

    return r


def get_gene(doc) -> str:
    """Get gene"""

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    
    try:

        namespace = doc["referenceEntity"]["databaseName"]
        id_ = doc["referenceEntity"]["identifier"]
        label = doc["referenceEntity"]["geneName"][0]
        label = clean_label(label)

    except Exception as e:
        id_ = ""
        namespace = ""
        label = ""
        errors_fh.write(f"{doc['dbId']}\tGene\t{str(e)}\n")
        log.exception(f"{doc['dbId']}\tGene\tError: {str(e)}\n")

    gene = Entity(namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid)
    g = Function(function="g", parameters=[gene], stid=stid, stid_version=stid_version, dbid=dbid)

    return g


def get_rna(doc) -> str:
    """Get RNA"""

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    
    try:
        prefix = doc["referenceEntity"]["databaseName"]
        _id = doc["referenceEntity"]["identifier"]
        label = doc["referenceEntity"]["geneName"][0]
        label = clean_label(label)

    except Exception as e:
        id_ = ""
        namespace = ""
        label = ""

        # errors_fh.write(f"{doc['dbId']}\RNA\t{str(e)}\n")
        log.exception(f"{doc['dbId']}\tRNA\tError: {str(e)}\n")

    rna = Entity(namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid)
    r = Function(function="r", parameters=[rna], stid=stid, stid_version=stid_version, dbid=dbid)

    return r


def get_polymer(doc) -> str:
    """Get Polymer"""

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    
    namespace = "REACTOME"
    stid = doc["stId"]
    label = doc["label"][0]
    label = clean_label(label)

    polymer = Entity(namespace, id=stid, label=label, stid=stid, stid_version=stid_version, dbid=dbid)
    p = Function(function="FNTBDPolymer", parameters=[polymer], stid=stid, stid_version=stid_version, dbid=dbid)
    return p


def get_other_entity(doc) -> str:
    """Get OtherEntity"""


    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    
    namespace = "REACTOME"
    stid = doc["stId"]
    label = doc["name"][0]
    label = clean_label(label)

    polymer = Entity(namespace, id=stid, label=label, stid=stid, stid_version=stid_version, dbid=dbid)
    p = Function(function="FNTBDOtherEntity", parameters=[polymer], stid=stid, stid_version=stid_version, dbid=dbid)
    return p


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


def process_candidate_set(candidateset):

    # print("Processing CandidateSet")
    dbid = candidateset.dbid
    stid = candidateset.stid
    label = clean_label(candidateset.label)
    fn = get_function_type(dbid, candidateset.members)
    sets[dbid] = []
    for member in candidateset.members:
        sets[dbid].append(
            {"subject": f"{fn}(REACTOME:{stid}!{label})", "relation": "isA", "object": member,}
        )

    # debug(sets[dbid])
    return f"{fn}(REACTOME:{stid}!{label})"


def format_complex(complex, complex_list: list, level: int = 0):
    """Format complex"""

    global sets

    # print("Start Complex", complex)

    normal_members = []
    set_members = []

    for member in complex.members:
        if isinstance(member, Complex):
            # print(f"{' ' * level* 2}Complex member {member.dbid}")
            level += 1
            complex_list.extend(format_complex(member, complex_list, level))
        elif isinstance(member, EntitySet):
            if member.type == "CandidateSet":
                candidateset = process_candidate_set(member)
                # print(f"{' ' * level * 2}CandidateSet {member.dbid}", candidateset)
                normal_members.append(candidateset)

            else:
                # print(f"{' ' * level * 2}DefinedSet {member.dbid}")
                set_members.append(member.members)
        else:
            # print(f"{' ' * level * 2}Member")
            normal_members.append(member)

    temp_list = []
    enumerated_set_members = list(itertools.product(*set_members))
    if enumerated_set_members:
        for esm in enumerated_set_members:
            temp_list.append(list(esm) + normal_members)
    else:
        temp_list = normal_members

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


def create_assertions(doc):

    comps = get_reaction_components(doc)
    assertions = []

    # print("DumpVar:\n", json.dumps(comps, cls=CustomEncoder, indent=4))
    
    catalysts = comps.get("catalysts", [])
    regulators = comps.get("regulators", [])
    inputs = comps.get("inputs", [])
    outputs = comps.get("outputs", [])

    if len(regulators) > 0 and len(catalysts) > 0:
        for regulator in regulators:
            for catalyst in catalysts:
                regulator_str = regulator.regulator.to_bel()
                catalyst_str = catalyst.catalyst.to_bel()
                assertion = f"{regulator_str} {regulator.relation} {catalyst_str}"
                assertions.append(assertion)

    elif len(regulators) > 0 and len(catalysts) == 0:
        with open(f"checks/missing_catalyst_{comps['dbid']}.json", 'w') as f:
            json.dump(comps, f, indent = 4)

    # Process input and outputs into a BEL construct
    target = ""
    check_inputs, check_outputs = [], []
    for input in inputs:
        if input.function == "a":
            continue
        check_inputs.append(input)
    for output in outputs:
        if output.function == "a":
            continue
        check_outputs.append(output)
    
    if len(check_inputs) == 1 and len(check_outputs) == 1:
        if check_inputs[0].loc != check_outputs[0].loc:  # Is this a translocation?
            target = f"tloc({check_inputs[0].to_bel()}, {check_outputs[0].to_bel()})"
        elif not check_inputs[0].mods and check_outputs[0].mods:  # Is this a modification?
            target = check_outputs[0].to_bel()
        else:  # default to a rxn()
            reactants = ", ".join([input.to_bel() for input in inputs])
            products = ", ".join([output.to_bel() for output in outputs])

            target = f"rxn(reactants({reactants}), products({products}))"

    if catalysts:
        for catalyst in catalysts:
            assertions.append(f"{catalyst.catalyst.to_bel()} {catalyst.relation} {target}")
    else:
        assertions.append(target)

    return assertions


def create_nanopub_from_reaction(reaction_id):

    doc = db.get_entity(reaction_id)

    if doc["schemaClass"] == "FailedReaction":
        # TODO - handle Failed Reaction - e.g. 3229118 (http://localhost:18529/_db/reactome_db/_admin/aardvark/index.html#collection/reactome/3229118)
        return {}

    src_url = f"https://reactome.org/content/detail/{doc['stIdVersion']}"
    metadata = {
        "source_url": src_url,
        "source": "Reactome",
        "license": "CC0",
        "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
        "gd_updateTS": now,
    }

    log.info("Reaction", dbid=doc["dbId"])

    metadata = get_creator(metadata, doc)

    nanopub = {"assertions": [], "annotations": []}
    annotations = []
    nanopub["id"] = f"Reactome_{doc['stIdVersion']}"

    nanopub["evidence"] = doc["summation"][0]["text"]
    nanopub["citation"] = {"uri": src_url}

    if "species" in doc:
        species_id = doc["species"][0]
        if isinstance(species_id, dict):
            species_id = species_id["dbId"]

        annotations.append(get_species(species_id))  # TODO - can there be multiple species?

    if "disease" in doc:
        for disease_obj in doc["disease"]:
            annotations.append(get_disease(disease_obj["dbId"]))

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
            nanopub = create_nanopub_from_reaction(reaction_id)

            print("Nanopub:\n", json.dumps(nanopub, indent=4))

            if not nanopub:
                continue

            # reactions.append(nanopub)
            f.write(f"{json.dumps(nanopub, cls=CustomEncoder)}\n")

    except Exception as e:
        log.exception("Converting reactions", error=e)
    finally:
        f.close()
