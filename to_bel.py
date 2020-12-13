# Standard Library
import copy
import datetime
import glob
import itertools
import json
import re
from collections.abc import Iterable
from typing import List, Optional

# Third Party Imports
import boltons.iterutils
import cachetools
import cachetools.func
import xxhash
from bel.core.utils import namespace_quoting
from bel.lang.ast import AssertionStr, BELAst
from devtools import debug, pprint

# Local Imports
import db
import settings
from db import collect_all_reactions
from logging_setup import log
from utils import clean_label, ts_convert

# Use new BEL entity format - <NS>:<ID>!<LABEL>, e.g. CHEBI:15999!water
NEW_BEL_ENTITY = True

db_objects = db.get_db()
reactome_db = db_objects["reactome_db"]
reactome_db_name = db_objects["reactome_db_name"]
reactome_coll = db_objects["reactome_coll"]
reactome_coll_name = db_objects["reactome_coll_name"]
reactome_aql = db_objects["reactome_aql"]

unmatched_fh = open("unmatched.tsv", "w")
errors_fh = open("errors.tsv", "w")

entity_sets = {}  # named families to write into BEL statemens (isA relations)
complexes = {}  # named complexes to write into BEL statements (hasComponent relations)

now = f"{datetime.datetime.utcnow().isoformat(timespec='milliseconds')}Z"


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
            entity_id = namespace_quoting(self.id)
            return f"{self.namespace}:{entity_id}"
        else:
            entity_id = namespace_quoting(self.id)
            entity_label = namespace_quoting(self.label)

            if NEW_BEL_ENTITY:
                entity = f"{self.namespace}:{entity_id}!{entity_label}"
            else:
                entity = f"{self.namespace}:{entity_id}"

            return entity

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.to_bel()


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

    def to_bel(self, strip_locations=False, add_location=True, complex_member=False):

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
                if not mod:
                    continue
                modifications.append(mod.to_bel())

            modifications = f', {", ".join(sorted(modifications))}'
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

    def __hash__(self):
        return int(xxhash.xxh32(str(self.__repr__())).intdigest())


class PMod:
    def __init__(
        self,
        modification,
        dbid,
        residue: Optional[str] = None,
        coordinate: Optional[int] = None,
        name: str = None,
    ):
        self.modification = modification
        self.residue = residue
        self.coordinate = coordinate
        self.dbid = dbid
        self.name = name

        if self.name:
            if NEW_BEL_ENTITY:
                self.modification = f"{self.modification}!{namespace_quoting(self.name)}"
            else:
                self.modification = f"{self.modification}"

    def to_bel(self, strip_locations=False):

        if self.modification and self.residue and self.coordinate:
            return f"pmod({self.modification}, {self.residue}, {self.coordinate})"
        # elif self.modification and self.coordinate:
        #     return f"pmod({self.modification}, , {self.coordinate})"
        # elif self.modification and self.residue:
        #     return f"pmod({self.modification}, {self.residue})"
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


class Variance:
    def __init__(
        self,
        variance,
        dbid,
    ):
        self.variance = variance
        self.dbid = dbid

    def to_bel(self, strip_locations=False):

        return f'var("p.{self.variance}")'

    def mod_only_to_bel(self):
        """Used for EntitySets"""
        return self.to_bel(self)

    def to_json(self):
        r = {
            "var": {
                "variance": self.variance,
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
        self.mods: List[PMod] = []
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

    def get_location(self):

        locations = set()
        for member in self.members:
            if member.loc in locations:
                # log.debug("Complex location - skipping", location=member.loc)
                continue
            if member.type == "Function" and member.loc:
                locations.add(member.loc)

        # log.debug("Locations", locations=locations)

        if len(locations) == 1:
            self.loc = list(locations)[0]
        elif len(locations) == 0:
            pass
        else:
            # TODO - figure out how to handle multiple locations
            #   [warning  ] Complex has mis-matched locations complex=190418 file=./to_bel.py function=get_location line=310 locations={'GO:0005576!"extracellular region"', 'GO:0005886!"plasma membrane"'}
            #   [warning  ] Complex has mis-matched locations complex=909708 file=./to_bel.py function=get_location line=310 locations={'GO:0005886!"plasma membrane"', 'GO:0005829!cytosol'}

            # log.warning("Complex has mis-matched locations", locations=locations, complex=self.dbid)
            self.loc = ""

    def get_members(self):
        """Get Complex members to flatten out the complex members"""

        members = set()
        for member in self.members:
            if member in members:
                log.debug("Complex skipping member", member=member)
                continue
            if isinstance(member, Complex):
                for c in member.get_members():
                    members.add(c)
            else:
                members.add(member)

        return list(members)

    def to_bel(self, strip_locations=True, add_location=True, complex_member: bool = False):
        strip_locations = True  # override this variable

        if add_location:
            self.get_location()
        else:
            self.loc = ""

        members_strings = []
        for member in sorted(self.get_members()):
            members_strings.append(
                member.to_bel(strip_locations=strip_locations, complex_member=True)
            )

        log.debug("Complex members", complex_member=complex_member, members_strings=members_strings)

        if complex_member:
            return ", ".join(members_strings)
        else:
            if self.loc:
                return f"complex({', '.join(members_strings)}, loc({self.loc}))"
            else:
                return f"complex({', '.join(members_strings)})"

    def __str__(self):
        return self.to_bel()

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return self.__repr__() < other.__repr__()

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    def __hash__(self):
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
        namespace: str = "reactome",
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
        self.function = "MISSING"
        self.loc = loc
        self.mods: List[PMod] = []

        self.get_location()  # Update self.loc
        self.get_function()  # Update self.function
        self.get_modifications()  # Update self.mods

        log.debug("EntitySet function", function=self.function, dbid=self.dbid)

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
            log.debug
            if isinstance(member, (Function, Complex)):
                functions.add(member.function)
            else:
                log.debug(f"EntitySet {self.dbid} has non-function member: {str(member)}")

        if len(functions) == 1:
            self.function = list(functions)[0]
        elif "complex" in functions:
            self.function = "complex"
        else:
            log.debug(f"EntitySet {self.dbid} has multiple or no functions {str(list(functions))}")

    def get_location(self):
        locations = set()
        for member in self.members:
            if isinstance(member, Function):
                locations.add(member.loc)
            else:
                log.debug(f"EntitySet {self.dbid} has no location: {str(member)}")

        if len(locations) == 1:
            self.loc = list(locations)[0]
        else:
            log.debug(f"EntitySet {self.dbid} has multiple or no locations {str(list(locations))}")

    def get_modifications(self):
        mods = set()
        for member in self.members:
            if isinstance(member, Function):
                for mod in member.mods:
                    mods.add(mod)

        if len(mods) > 0:
            self.mods = list(mods)
        else:
            self.mods = []

        self.mods = [mod for mod in self.mods if mod is not None]

    def to_bel(self, strip_locations=False, add_location=True, complex_member=False):

        location = ""
        if not strip_locations and self.loc:
            location = f", loc({self.loc})"

        modifications = ""

        if self.mods:
            modifications = list(set([mod.mod_only_to_bel() for mod in self.mods]))
            modifications = ", ".join(modifications)
            modifications = f", {modifications}"

        if complex_member and self.function == "complex":
            if NEW_BEL_ENTITY:
                entity = f"reactome:{self.stid_version}!{namespace_quoting(self.label)}{location}"
            else:
                entity = f"reactome:{self.stid_version}{location}"

            return entity

        else:
            if NEW_BEL_ENTITY:
                entity = f"{self.function}(reactome:{self.stid_version}!{namespace_quoting(self.label)}{modifications}{location})"
            else:
                entity = f"{self.function}(reactome:{self.stid_version}{modifications}{location})"
            return entity

    def list_members(self):
        members_list = []
        for member in self.members:
            if isinstance(member, str):
                members_list.append(member)
            else:
                members_list.append(str(member))

        return f'{self.type}({", ".join(members_list)}, loc({self.loc}))'

    def __str__(self):
        if NEW_BEL_ENTITY:
            return f"reactome:{self.stid_version}!{namespace_quoting(self.label)} - {self.type} - dbid: {self.dbid}"
        else:
            return f"reactome:{self.stid_version} - {self.type} - dbid: {self.dbid}"

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

    def __hash__(self):
        return int(xxhash.xxh32(str(self.__repr__())).intdigest())


class Regulator:
    def __init__(
        self,
        regulator,
        relation,
        dbid,
        stid,
        stid_version: str,
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
        self,
        catalyst,
        relation,
        dbid,
        stid,
        stid_version: str,
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


@cachetools.func.lru_cache(maxsize=1024, typed=False)
def get_creator(dbid):

    doc = db.get_entity(dbid)
    metadata = {}

    if "authored" in doc:
        if isinstance(doc["authored"], list):
            authored_id = doc["authored"][0]
        else:
            authored_id = doc["authored"]

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


@cachetools.func.lru_cache(maxsize=1024, typed=False)
def get_species(species_id):
    doc = db.get_entity(species_id)
    return {"type": "Species", "id": f"TAX:{doc['taxId']}", "label": doc["displayName"]}


@cachetools.func.lru_cache(maxsize=1024, typed=False)
def get_disease(dbid):

    doc = db.get_entity(dbid)

    # import json
    # print('DumpVar:\n', json.dumps(doc, indent=4))

    if doc["databaseName"] == "DOID":
        key = f"DO:{doc['identifier']}"
    else:
        key = f"{doc['databaseName']}:{doc['identifier']}"

    return {
        "type": "Disease",
        "id": key,
        "label": doc["displayName"],
    }


@cachetools.func.lru_cache(maxsize=1024, typed=False)
def process_mod(dbid):
    doc = db.get_entity(dbid)

    # process pmod()
    if dbid == 917934:
        return PMod(
            modification="Ac",
            residue="A",
            coordinate=2,
            dbid=dbid,
            name="N-acetyl-L-alanine",
        )

    elif doc.get("className", "") == "FragmentDeletionModification":
        start = doc.get("startPositionInReferenceSequence", None)
        end = doc.get("endPositionInReferenceSequence", None)

        if start and end:
            variance = f"p.{start}_{end}del"
            return Variance(variance, dbid)

    elif "psiMod" in doc and doc.get("className", "") == "ReplacedResidue":
        variance = convert_replace_residue_display_name(doc.get("displayName", ""))
        return Variance(variance, dbid)

    elif "psiMod" in doc:

        coordinate = doc.get("coordinate", "")

        if isinstance(doc["psiMod"], dict):
            psimods = [doc["psiMod"]]
        else:
            psimods = doc["psiMod"]

        for psimod in psimods:

            modification = ""
            residue = ""
            name = ""
            if psimod["identifier"] == "00046":
                modification = "Ph"
                residue = "Ser"

            elif psimod["identifier"] == "00047":
                modification = "Ph"
                residue = "Thr"

            else:
                modification = f'PSIMOD:{psimod["identifier"]}'
                if "name" in psimod:
                    name = psimod["name"][0]

            return PMod(
                modification=modification,
                residue=residue,
                coordinate=coordinate,
                dbid=dbid,
                name=name,
            )

    else:
        log.error(
            f"Unable to process protein modification: {dbid}  classname: {doc['className']}  displayName: {doc['displayName']}"
        )
        return PMod("Missing", dbid)


def convert_replace_residue_display_name(display_name):

    # print("Display name", display_name)
    pattern = f"({settings.amino_acids_regex}).*?\s+(\d+)\s+.*?({settings.amino_acids_regex})"
    # print("Pattern", pattern)
    match = re.search(
        pattern,
        display_name,
        flags=re.IGNORECASE,
    )

    # print("Groups", match.groups(), match.group(1))

    start_residue = match.group(1).lower()
    position = match.group(2).lower()
    end_residue = match.group(3).lower()

    variance = f"{settings.amino_acids[start_residue].title()}{position}{settings.amino_acids[end_residue].title()}"

    return variance


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_protein(dbid) -> Entity:

    doc = db.get_entity(dbid)

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

        # log.info("Get Label for protein", doc=doc)

        name = ""
        if "name" in doc["referenceEntity"] and len(doc["referenceEntity"]["name"]) > 0:
            name = doc["referenceEntity"]["name"][0]
        elif "geneName" in doc["referenceEntity"] and len(doc["referenceEntity"]["geneName"]) > 0:
            name = doc["referenceEntity"]["geneName"][0]

        label = clean_label(name)

        id_ = ""
        if (
            "identifier" in doc["referenceEntity"]
            and isinstance(doc["referenceEntity"]["identifier"], list)
            and len(doc["referenceEntity"]["identifier"]) > 0
        ):
            id_ = doc["referenceEntity"]["identifier"][0]
        elif "identifier" in doc["referenceEntity"]:
            id_ = str(doc["referenceEntity"]["identifier"])

        loc = ""
        if doc.get("compartment", False):  # TODO - can there be multiple compartments?
            if doc["compartment"][0]["databaseName"] == "GO":
                loc_namespace = "GO"
            else:
                log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
                loc_namespace = "TBD"

            loc_id = doc["compartment"][0]["accession"]
            loc_label = clean_label(doc["compartment"][0]["displayName"])
            if NEW_BEL_ENTITY:
                loc = f"{loc_namespace}:{loc_id}!{loc_label}"
            else:
                loc = f"{loc_namespace}:{loc_id}"

        if "hasModifiedResidue" in doc:
            for modification in doc["hasModifiedResidue"]:
                if isinstance(modification, int):
                    modification = db.get_entity(modification)
                mods.append(process_mod(modification["dbId"]))

    except Exception as e:
        errors_fh.write(f"{doc['dbId']}\tProtein\t{str(e)}\n")
        log.exception(f"{doc['dbId']}\tProtein\tError: {str(e)}\n")
        namespace = "reactome"
        id_ = dbid
        label = ""

    protein = Entity(
        namespace=namespace,
        id=id_,
        label=label,
        stid=stid,
        stid_version=stid_version,
        dbid=dbid,
        loc=loc,
    )
    p = Function(
        function="p",
        parameters=[protein],
        mods=mods,
        stid=stid,
        stid_version=stid_version,
        dbid=dbid,
        loc=loc,
    )

    return p


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_catalyst(catalyst_dbid) -> str:
    """Get catalyst"""

    doc = db.get_entity(catalyst_dbid)

    try:
        if isinstance(doc["physicalEntity"], int):
            dbid = doc["physicalEntity"]
        elif isinstance(doc["physicalEntity"], dict):
            dbid = doc["physicalEntity"]["dbId"]
        else:
            log.error("Cannot determine dbId for 'physicalEntity' in reaction_id: {dbid}")

        doc = db.get_entity(dbid)
        stid = ["stId"]
        stid_version = ["stIdVersion"]

        catalyst = process_component(dbid)
        if catalyst and catalyst.function != "a":
            catalyst = Function(
                function="act",
                parameters=[catalyst],
                stid=stid,
                stid_version=stid_version,
                dbid=dbid,
            )

        # relation
        relation = "directlyIncreases"

        r = Catalyst(catalyst, relation=relation, dbid=dbid, stid=stid, stid_version=stid_version)

    except Exception as e:
        errors_fh.write(f"Catalyst dbid: {dbid}\tCatalyst\t{str(e)}\n")
        log.error(f"No dbid for physical entity in {catalyst_dbid}\tCatalyst\tError: {str(e)}\n")
        catalyst = f"Cannot process Catalyst for {catalyst_dbid}"
        r = Catalyst(catalyst, relation="", dbid="", stid="", stid_version="")

    return r


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_gt(dbid) -> str:
    """Get className Genes and Transcripts - seems like it's used for proteins

    Example: http://thor:9529/_db/_system/_admin/aardvark/index.html#collection/reactome/381096

    """

    doc = db.get_entity(dbid)

    log.debug(f"GT className: {doc['dbId']}")

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    loc = ""

    if "referenceEntity" in doc:
        label = clean_label(doc["referenceEntity"]["name"][0])
        if doc["referenceEntity"]["databaseName"] == "UniProt":
            id_: str = doc["referenceEntity"]["identifier"]
            namespace = "SP"
    else:
        namespace = "TBD"
        label = clean_label(doc["name"][0])
        id_ = label

    loc = None
    if doc.get("compartment", False):  # TODO - can there be multiple compartments?
        if doc["compartment"][0]["databaseName"] == "GO":
            loc_prefix = "GO"
        else:
            log.info(f'Missing location databaseName {doc["compartment"][0]["databaseName"]}')
            loc_prefix = "TBD"

        loc_id = doc["compartment"][0]["accession"]
        loc_name = clean_label(doc["compartment"][0]["displayName"])
        if NEW_BEL_ENTITY:
            loc = f"{loc_prefix}:{loc_id}!{loc_name}"
        else:
            loc = f"{loc_prefix}:{loc_id}"

    protein = Entity(
        namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid, loc=loc
    )
    p = Function("p", [protein], stid, stid_version, dbid, loc=loc)

    return p


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_compound(dbid: str) -> str:

    doc = db.get_entity(dbid)

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
        if NEW_BEL_ENTITY:
            loc = f"{loc_prefix}:{loc_id}!{loc_name}"
        else:
            loc = f"{loc_prefix}:{loc_id}"

    compound = Entity(
        namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid, loc=loc
    )
    c = Function("a", [compound], stid, stid_version, dbid, loc=loc)

    return c


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_complex(dbid: str) -> str:
    """Get complex"""

    doc = db.get_entity(dbid)

    components = set()
    for component in doc["hasComponent"]:
        if isinstance(component, int):
            dbid = component
        else:
            dbid = component["dbId"]

        r = process_component(dbid)
        if r:
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

        if NEW_BEL_ENTITY:
            loc = f"{loc_prefix}:{loc_id}!{loc_name}"
        else:
            loc = f"{loc_prefix}:{loc_id}"

    c = Complex(
        list(components), label=label, dbid=dbid, stid=stid, stid_version=stid_version, loc=loc
    )

    # Adding
    complexes[dbid] = c

    return c


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_set(dbid: str) -> dict:

    doc = db.get_entity(dbid)
    # json.dump(doc, sets_fh, indent=4)
    # log.info("Set", doc=doc)

    label = doc["name"][0]
    members = set()

    if "hasMember" in doc:
        for member in doc["hasMember"]:
            r = process_component(member["dbId"])
            if r:
                members.add(r)

    if "hasCandidate" in doc:
        # label += " [CANDIDATE SET]"
        for member in doc["hasCandidate"]:
            r = process_component(member["dbId"])
            if r:
                members.add(r)

    label = clean_label(label)

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

        if NEW_BEL_ENTITY:
            loc = f"{loc_prefix}:{loc_id}!{loc_name}"
        else:
            loc = f"{loc_prefix}:{loc_id}"

    s = EntitySet(
        members=list(members),
        namespace="reactome",
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


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_regulator(regulator_dbid: str) -> str:
    """Get regulation"""

    regulator_doc = db.get_entity(regulator_dbid)

    # print("DumpVar:\n", json.dumps(doc, indent=4))

    if "regulator" not in regulator_doc:
        log.error("Bad regulator doc - missing regulator attribute", dbid=dbid)
        return Regulator("", relation="", dbid="", stid="", stid_version="")

    try:
        if isinstance(regulator_doc["regulator"], int):
            dbid = regulator_doc["regulator"]
        elif isinstance(regulator_doc["regulator"], dict):
            dbid = regulator_doc["regulator"]["dbId"]
        else:
            log.error(f"No dbid in regulator doc: {regulator_dbid}")

        doc = db.get_entity(dbid)

        stid = ["stId"]
        stid_version = doc["stIdVersion"]

        regulator = process_component(dbid)
        if not regulator:
            return None

        # Wrap as act() function
        if regulator.function != "a":
            regulator = Function(
                function="act",
                parameters=[regulator],
                stid=stid,
                stid_version=stid_version,
                dbid=dbid,
            )

        # relation
        if regulator_doc["className"] in [
            "PositiveGeneExpressionRegulation",
            "PositiveRegulation",
        ]:
            relation = "directlyIncreases"

        elif regulator_doc["className"] in [
            "NegativeGeneExpressionRegulation",
            "NegativeRegulation",
        ]:
            relation = "directlyDecreases"
        else:
            log.info(f"Unknown regulation relationship: {regulator_dbid}")
            relation = "regulates"

        r = Regulator(regulator, relation=relation, dbid=dbid, stid=stid, stid_version=stid_version)

        return r

    except Exception as e:
        errors_fh.write(f"{regulator_dbid}\tRegulator\t{str(e)}\n")
        log.exception(f"{regulator_dbid}\tRegulator\tError: {str(e)}\n")

        return None


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_gene(dbid) -> str:
    """Get gene"""

    doc = db.get_entity(dbid)

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


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_rna(dbid) -> str:
    """Get RNA"""

    doc = db.get_entity(dbid)
    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]
    label = ""
    namespace = ""
    id_ = ""

    try:
        namespace = doc["referenceEntity"]["databaseName"]
        id_ = doc["referenceEntity"]["identifier"]
        if "genename" in doc["referenceEntity"] and len(doc["referenceEntity"]["geneName"]) >= 0:
            label = doc["referenceEntity"]["geneName"][0]
        label = clean_label(label)

    except Exception as e:
        log.exception(f"{doc['dbId']}\tRNA\tError: {str(e)}\n")

    rna = Entity(namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid)
    r = Function(function="r", parameters=[rna], stid=stid, stid_version=stid_version, dbid=dbid)

    return r


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_polymer(dbid) -> str:
    """Get Polymer"""

    doc = db.get_entity(dbid)

    stid = doc["stId"]
    stid_version = doc["stIdVersion"]
    dbid = doc["dbId"]

    namespace = "reactome"
    stid = doc["stId"]
    if "label" in doc:
        label = doc["label"][0]
    elif "name" in doc:
        label = doc["name"][0]
    else:
        label = ""

    label = clean_label(label)

    polymer = Entity(
        namespace, id=stid, label=label, stid=stid, stid_version=stid_version, dbid=dbid
    )
    p = Function(
        function="FNTBDPolymer",
        parameters=[polymer],
        stid=stid,
        stid_version=stid_version,
        dbid=dbid,
    )
    return p


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_other_entity(dbid) -> str:
    """Get OtherEntity"""

    doc = db.get_entity(dbid)
    stid = doc["stId"]
    stid_version = doc["stIdVersion"]

    if doc["name"][0] == "Photon":
        e = Entity(
            namespace="CHEBI",
            id="30212",
            label="photon",
            stid=stid,
            dbid=dbid,
            stid_version=stid_version,
        )

        a = Function("a", [e], stid, stid_version, dbid)
        return a

    log.warning(f"FNTBDOtherEntity dbid={dbid}  name={doc['name']}")
    return None

    # doc = db.get_entity(dbid)

    # stid = doc["stId"]
    # stid_version = doc["stIdVersion"]
    # dbid = doc["dbId"]

    # namespace = "reactome"
    # stid = doc["stId"]
    # label = doc["name"][0]
    # label = clean_label(label)

    # polymer = Entity(
    #     namespace, id=stid, label=label, stid=stid, stid_version=stid_version, dbid=dbid
    # )
    # p = Function(
    #     function="FNTBDOtherEntity",
    #     parameters=[polymer],
    #     stid=stid,
    #     stid_version=stid_version,
    #     dbid=dbid,
    # )
    # log.warning("FNTBDOtherEntity", dbid=dbid)
    # return p


@cachetools.func.lru_cache(maxsize=5000, typed=False)
def get_requirement(dbid):

    doc = db.get_entity(dbid)

    if isinstance(doc["regulator"], int):
        dbid = doc["regulator"]
    else:
        dbid = doc["regulator"]["dbId"]

    doc = db.get_entity(dbid)

    if doc["className"] == "Complex":
        regulator = get_complex(dbid)
    elif doc["className"] == "Chemical Compound":
        regulator = get_compound(dbid)
    elif doc["className"] == "Protein":
        regulator = get_protein(dbid)
    elif doc["className"] == "Set":
        regulator = get_set(dbid)
    else:
        log.error(f"Unknown class for Requirement regulator: {doc['className']}")
        regulator = None

    if regulator:
        stid = doc["stId"]
        stid_version = doc["stIdVersion"]

        return Regulator(
            regulator, relation="increases", dbid=dbid, stid=stid, stid_version=stid_version
        )


def process_component(dbid: str) -> dict:
    """Process reaction component: input, output, catalyst or regulator"""

    doc = db.get_entity(dbid)

    if doc["className"] == "Protein":
        return get_protein(dbid)
    elif doc["className"] == "Complex":
        return get_complex(dbid)
    elif doc["className"] == "Requirement":
        return get_requirement(dbid)
    elif doc["className"] in ["Chemical Compound", "ChemicalDrug", "ProteinDrug"]:
        return get_compound(dbid)
    elif doc["className"] == "Set":
        return get_set(dbid)
    elif doc["className"] == "DNA Sequence":
        return get_gene(dbid)
    elif doc["className"] == "RNA Sequence":
        return get_rna(dbid)
    elif doc["className"] == "Genes and Transcripts":
        return get_gt(dbid)
    elif doc["className"] == "Polymer":
        return get_polymer(dbid)
    elif doc["className"] == "OtherEntity":
        return get_other_entity(dbid)
    elif doc["className"] == "CatalystActivity":
        return get_catalyst(dbid)
    elif doc["className"] in [
        "PositiveGeneExpressionRegulation",
        "NegativeGeneExpressionRegulation",
        "NegativeRegulation",
        "PositiveRegulation",
    ]:
        return get_regulator(dbid)
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

    inputs = []
    outputs = []
    catalysts = []
    regulators = []

    reaction_type = doc.get("category", "")

    # inputs
    for input in doc.get("input", []):
        if isinstance(input, int):
            dbid = input
        else:
            dbid = input["dbId"]
        r = process_component(dbid)
        if r:
            inputs.append(r)

    for output in doc.get("output", []):
        if isinstance(output, int):
            dbid = output
        else:
            dbid = output["dbId"]
        r = process_component(dbid)
        if r:
            outputs.append(r)

    for regulator in doc.get("regulatedBy", []):
        if isinstance(regulator, int):
            dbid = regulator
        else:
            dbid = regulator["dbId"]
        r = process_component(dbid)
        if r:
            regulators.append(r)

    for catalyst in doc.get("catalystActivity", []):
        if isinstance(catalyst, int):
            dbid = catalyst
        else:
            dbid = catalyst["dbId"]

        r = process_component(dbid)
        if r:
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


def create_assertions(doc):

    comps = get_reaction_components(doc)
    assertions = []

    # print("DumpVar:\n", json.dumps(comps, cls=CustomEncoder, indent=4))

    catalysts = comps.get("catalysts", [])
    regulators = comps.get("regulators", [])
    inputs = comps.get("inputs", [])
    outputs = comps.get("outputs", [])

    # Add regulator -> catalyst Assertions
    #    These are created before the reaction assertions
    if len(regulators) > 0 and len(catalysts) > 0:
        for regulator in regulators:
            if not regulator:
                continue
            for catalyst in catalysts:
                if not catalyst:
                    continue
                regulator_str = regulator.regulator.to_bel()
                catalyst_str = catalyst.catalyst.to_bel()
                assertions.append(
                    {
                        "subject": regulator_str,
                        "relation": regulator.relation,
                        "object": catalyst_str,
                    }
                )

    # Create Assertion Objects (target)
    # Check for protein expression, protein modification or translocation assertions

    target = ""  # the object of the assertion

    check_inputs, check_outputs = [], []
    for input in inputs:
        if input.function == "a":  # ignore ATP/GTP/etc for modifications/translocations
            continue
        check_inputs.append(input)

    for output in outputs:
        if output.function == "a":  # ignore ATP/GTP/etc for modifications/translocations
            continue
        check_outputs.append(output)

    log.debug("Check input/outputs", check_inputs=check_inputs, check_outputs=check_outputs)

    if len(check_inputs) == 1 and len(check_outputs) == 1:
        check_input = check_inputs[0]
        check_output = check_outputs[0]

        # Gene expression -> p()
        if check_input.function == "g" and check_output.function == "p":
            target = check_output.to_bel()

        # Protein modification
        elif not check_input.mods and check_output.mods:
            target = check_output.to_bel()

        # TODO - check for modified complex - e.g. Reaction dbid = 1981128

        # Translocation expression
        elif check_input.loc and check_output.loc and check_input.loc != check_output.loc:
            from_loc = check_input.loc
            to_loc = check_output.loc
            target = f"tloc({check_input.to_bel(add_location=False, strip_locations=True)}, fromLoc({from_loc}), toLoc({to_loc}))"

    # default to a rxn() object
    if not target:
        reactants = ", ".join(sorted([input.to_bel() for input in inputs]))
        products = ", ".join(sorted([output.to_bel() for output in outputs]))

        target = f"rxn(reactants({reactants}), products({products}))"

    # If only catalysts - <catalyst> <catalyst.relation> <target>
    if catalysts:
        for catalyst in catalysts:
            assertions.append(
                {
                    "subject": catalyst.catalyst.to_bel(),
                    "relation": catalyst.relation,
                    "object": target,
                }
            )

    # If only regulators - <regulator> <regulator.relation> <target>
    elif len(catalysts) == 0 and len(regulators) > 0:
        for regulator in regulators:
            assertions.append(
                {
                    "subject": regulator.regulator.to_bel(),
                    "relation": regulator.relation,
                    "object": target,
                }
            )

    # No subject - only reaction
    elif target.startswith("rxn") or target.startswith("tloc"):
        assertions.append({"subject": target})

    for assertion in assertions:
        assertion_obj = AssertionStr(
            subject=assertion["subject"],
            relation=assertion.get("relation", ""),
            object=assertion.get("object", ""),
        )
        ast = BELAst(assertion=assertion_obj)
        ast.validate()
        for error in ast.errors:
            if error.severity == "error":
                log.error(f"Assertion: {assertion}  Error: {str(error)}")

    return assertions


def create_nanopub_from_reaction(reaction_id):

    doc = db.get_entity(reaction_id)

    if doc["schemaClass"] == "FailedReaction":
        log.warning(f"FailedReaction - reaction_id: {reaction_id}")
        # TODO - handle Failed Reaction - e.g. 3229118 (http://localhost:18529/_db/reactome_db/_admin/aardvark/index.html#collection/reactome/3229118)
        return {"errors": ["Failed Reaction"]}

    src_url = f"https://reactome.org/content/detail/{doc['stIdVersion']}"
    metadata = {
        "source_url": src_url,
        "source": "Reactome",
        "license": "CC0",
        "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
        "gd_updateTS": now,
        "collections": ["Reactome"],
        "version": f"Reactome_{doc['stIdVersion']}",
    }

    metadata.update(get_creator(reaction_id))

    nanopub = {"type": {"name": "BEL", "version": "2.1.2"}, "assertions": [], "annotations": []}

    annotations = []
    nanopub["id"] = f"Reactome_{doc['stId']}"

    if "summation" in doc:
        nanopub["evidence"] = doc["summation"][0]["text"]
    nanopub["citation"] = {"uri": src_url}

    if "species" in doc:
        if isinstance(doc["species"], list):
            species_id = doc["species"][0]
        else:
            species_id = doc["species"]

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


def convert(reaction_ids: list = [], update_all: bool = False, limit: Optional[int] = None):
    """Convert reactome reactions to BEL

    If reaction_ids is empty - it will process all Reactions in database
    """

    test_flag = False
    if reaction_ids:
        test_flag = True

    if not reaction_ids:
        reaction_ids = collect_all_reactions()

    number_reactions = len(reaction_ids)

    f = open("reactome_nanopubs.jsonl", "w")

    converted_reaction_ids = []
    if not update_all and not test_flag:
        try:
            with open("converted_reaction_ids.json", "r") as fconverted:
                converted_reaction_ids = json.load(fconverted)
        except Exception as e:
            converted_reaction_ids = []

    if converted_reaction_ids:
        log.info(
            f"Skipping {len(converted_reaction_ids)} reactions already converted based on the converted_reaction_ids.json file"
        )

    # import pdb; pdb.set_trace()

    counter = 0
    for reaction_id in reaction_ids:
        if not update_all and converted_reaction_ids and reaction_id in converted_reaction_ids:
            continue

        counter += 1
        log.info(f"Starting Reaction {reaction_id}", counter=counter, total=number_reactions)

        try:
            nanopub = create_nanopub_from_reaction(reaction_id)

            if nanopub.get("errors", []) == ["Failed Reaction"]:
                pass  # Skip logging error for this

            elif len(nanopub.get("assertions", [])) == 0:
                log.error("No assertions", reaction=reaction_id)

            else:
                nanopub = {"nanopub": nanopub}
                f.write(f"{json.dumps(nanopub, cls=CustomEncoder)}\n")
            # print("Nanopub:\n", json.dumps(nanopub, indent=4))

            converted_reaction_ids.append(reaction_id)

            if counter % 100 == 0:
                log.info(f"Processed {counter} reactions")
                with open("converted_reaction_ids.json", "w") as fconverted:
                    json.dump(converted_reaction_ids, fconverted, indent=4)

        except Exception as e:
            nanopub = {}
            errors_fh.write(f"{reaction_id}\t{str(e)}\n")
            errors_fh.flush()
            log.error(f"Problem converting reaction {reaction_id}  error: {str(e)}")

        if limit and counter >= limit:
            break

    f.close()

    if not test_flag:
        with open("converted_reaction_ids.json", "w") as fconverted:
            json.dump(converted_reaction_ids, fconverted, indent=4)
