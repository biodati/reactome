from to_bel import create_nanopub_from_reaction
import json
from difflib import unified_diff
from pprint import pprint


def compare(expected, actual):
    """Compare two nanopubs"""

    # Remove elements that are always changing
    expected["metadata"].pop("gd_createTS", "")
    expected["metadata"].pop("gd_updateTS", "")
    actual["metadata"].pop("gd_createTS", "")
    actual["metadata"].pop("gd_updateTS", "")  

    expected["assertions"] = sorted(expected["assertions"])
    actual["assertions"] = sorted(actual["assertions"])

    expected_str = json.dumps(expected, sort_keys=True)
    actual_str = json.dumps(actual, sort_keys=True)

    if expected_str != actual_str:
        print("Expected", expected)
        print("Actual  ", actual)
        # pprint(list(unified_diff(expected_str, actual_str)))
        return False
    
    return True
    

def test_regulated_by():

    expected = {
        "assertions": [
            "a(CHEBI:16356!\"3',5'-cyclic GMP\", loc(GO:0005829!cytosol)) directlyDecreases activity(complex(a(CHEBI:29108!\"calcium(2+)\"), p(SP:P0DP23!CALM1), p(SP:Q15746!MYLK), loc(GO:0005829!cytosol))",
            "activity(complex(a(CHEBI:29108!\"calcium(2+)\"), p(SP:P0DP23!CALM1), loc(GO:0005829!cytosol)) directlyIncreases activity(complex(a(CHEBI:29108!\"calcium(2+)\"), p(SP:P0DP23!CALM1), p(SP:Q15746!MYLK), loc(GO:0005829!cytosol))",
            "activity(complex(a(CHEBI:29108!\"calcium(2+)\"), p(SP:P0DP23!CALM1), p(SP:Q15746!MYLK), loc(GO:0005829!cytosol)) directlyIncreases p(REACTOME:R-HSA-445770.1!\"Phosphorylated Smooth Muscle Myosin Light Chain\", pmod(Ph), loc(GO:0005829!cytosol))"
        ],
        "annotations": [
            {
                "label": "Homo sapiens",
                "type": "Species",
                "id": "TAX:9606",
                
            }
        ],
        "id": "Reactome_R-HSA-445813.2",
        "evidence": "The smooth muscle light chain kinase phosphorylates the smooth muscle light chains. This phosphorylation activates the myosin lights chains, effectively allowing contraction to begin.",
        "citation": {
            "uri": "https://reactome.org/content/detail/R-HSA-445813.2"
        },
        "metadata": {
            "source_url": "https://reactome.org/content/detail/R-HSA-445813.2",
            "source": "Reactome",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "gd_updateTS": "2020-02-28T16:12:38.245Z",
            "creator": "Gillespie, ME",
            "gd_createTS": "2009-11-18T21:36:06.000Z",
            "creator_orcid": "0000-0002-5766-1702"
        }
    }

    actual = create_nanopub_from_reaction("445813")
    assert actual

    result = compare(expected, actual)
    assert result


def test_disease():

    expected = {}

    actual = create_nanopub_from_reaction("2399988")
    assert actual

    result = compare(expected, actual)
    assert result


def test_catalyzed():

    expected = {}

    actual = create_nanopub_from_reaction("9027627")
    assert actual

    result = compare(expected, actual)
    assert result



def test_candidateset_subset():

    expected = {}

    actual = create_nanopub_from_reaction("5655336")
    assert actual

    result = compare(expected, actual)
    assert result


def test_nestedcomplexes_candidatesets():

    expected = {}

    actual = create_nanopub_from_reaction("8952044")
    assert actual

    result = compare(expected, actual)
    assert result


def test_nestedcomplexes_entitysets():

    expected = {}

    actual = create_nanopub_from_reaction("8942101")
    assert actual

    result = compare(expected, actual)
    assert result
