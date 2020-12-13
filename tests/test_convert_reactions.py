# Standard Library
import json

# Local Imports
from logging_setup import log
from to_bel import create_nanopub_from_reaction

# import structlog
# log = structlog.get_logger()


def compare(expected, actual):
    """Compare two nanopubs"""

    # Remove elements that are always changing
    if expected.get("metadata", False):
        expected["metadata"].pop("gd_createTS", "")
        expected["metadata"].pop("gd_updateTS", "")
    if actual.get("metadata", False):
        actual["metadata"].pop("gd_createTS", "")
        actual["metadata"].pop("gd_updateTS", "")

    expected["assertions"] = sorted(expected.get("assertions", []))
    actual["assertions"] = sorted(actual.get("assertions", []))

    expected_str = json.dumps(expected, sort_keys=True)
    actual_str = json.dumps(actual, sort_keys=True)

    expected_assertions_str = json.dumps(expected["assertions"], sort_keys=True)
    actual_assertions_str = json.dumps(actual["assertions"], sort_keys=True)

    # Dump out actual if no expected
    if not expected.get("id", False):
        print("\nActual formatted", json.dumps(actual, sort_keys=True, indent=4))

    if expected_str and expected_assertions_str != actual_assertions_str:
        print("\n\nASSERTIONS MISMATCH")
        print("\nExpected Assertions", expected["assertions"])
        print("\nActual    Assertions", actual["assertions"])

        return False

    elif expected_str != actual_str:
        print("\n\nNANOPUBS MISMATCH")
        print("\nExpected", expected_str)
        print("\nActual  ", actual_str)

        return False

    return True


def test_regulated_by():

    expected = {
        "assertions": [
            'activity(complex(a(CHEBI:29108!"calcium(2+)"), p(SP:P0DP23!CALM1), loc(GO:0005829!cytosol))) directlyIncreases activity(complex(a(CHEBI:29108!"calcium(2+)"), p(SP:P0DP23!CALM1), p(SP:Q15746!MYLK), loc(GO:0005829!cytosol)))',
            'a(CHEBI:16356!"3\',5\'-cyclic GMP", loc(GO:0005829!cytosol)) directlyDecreases activity(complex(a(CHEBI:29108!"calcium(2+)"), p(SP:P0DP23!CALM1), p(SP:Q15746!MYLK), loc(GO:0005829!cytosol)))',
            'activity(complex(a(CHEBI:29108!"calcium(2+)"), p(SP:P0DP23!CALM1), p(SP:Q15746!MYLK), loc(GO:0005829!cytosol))) directlyIncreases p(REACTOME:R-HSA-445770.1!"Phosphorylated Smooth Muscle Myosin Light Chain", pmod(Ph), loc(GO:0005829!cytosol))',
        ],
        "annotations": [{"type": "Species", "id": "TAX:9606", "label": "Homo sapiens"}],
        "id": "Reactome_R-HSA-445813.2",
        "evidence": "The smooth muscle light chain kinase phosphorylates the smooth muscle light chains. This phosphorylation activates the myosin lights chains, effectively allowing contraction to begin.",
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-445813.2"},
        "metadata": {
            "source_url": "https://reactome.org/content/detail/R-HSA-445813.2",
            "source": "Reactome",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "gd_updateTS": "2020-03-31T15:07:37.303Z",
            "creator": "Gillespie, ME",
            "gd_createTS": "2009-11-18T21:36:06.000Z",
            "creator_orcid": "0000-0002-5766-1702",
        },
    }

    actual = create_nanopub_from_reaction("445813")
    assert actual

    result = compare(expected, actual)
    assert result


def test_disease():

    expected = {
        "assertions": [
            'activity(p(SP:P31749!AKT1, pmod(PSIMOD:01636!"L-glutamic acid removal", , 17), pmod(Ph, Ser, 473), pmod(Ph, Thr, 308), loc(GO:0005654!nucleoplasm))) directlyIncreases p(SP:P22736!NR4A1, pmod(Ph, Ser, 351), loc(GO:0005654!nucleoplasm))'
        ],
        "annotations": [
            {"type": "Species", "id": "TAX:9606", "label": "Homo sapiens"},
            {"type": "Disease", "id": "DOID:162", "label": "cancer"},
        ],
        "id": "Reactome_R-HSA-2399988.1",
        "evidence": "AKT1 E17K gain-of-function mutant is expected to phosphorylate NR4A1, like the wild-type AKT (Pekarsky et al. 2001), but this has not been experimentally tested.",
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-2399988.1"},
        "metadata": {
            "source_url": "https://reactome.org/content/detail/R-HSA-2399988.1",
            "source": "Reactome",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "creator": "Orlic-Milacic, M",
            "creator_orcid": "0000-0002-3218-5631",
        },
    }
    actual = create_nanopub_from_reaction("2399988")
    assert actual

    result = compare(expected, actual)
    assert result


def test_catalyzed():

    expected = {
        "assertions": [
            'activity(complex(a(CHEBI:26355!"heme b"), p(SP:P35354!PTGS2, pmod(PSIMOD:00369!O-acetyl-L-serine, , 516)), loc(GO:0005789!"endoplasmic reticulum membrane"))) directlyIncreases rxn(reactants(a(CHEBI:15379!dioxygen, loc(GO:0005829!cytosol)), a(CHEBI:28125!"all-cis-docosa-4,7,10,13,16,19-hexaenoic acid", loc(GO:0005829!cytosol))), products(a(CHEBI:72637!"(4Z,7Z,10Z,13Z,15E,19Z)-17-hydroxydocosahexaenoic acid", loc(GO:0005829!cytosol))))'
        ],
        "annotations": [{"type": "Species", "id": "TAX:9606", "label": "Homo sapiens"}],
        "id": "Reactome_R-HSA-9027627.1",
        "evidence": "Normally, cyclooxygenases (COX) carry out stereospecific oxygenation of arachidonic acid to generate prostaglandins. When treated with aspirin (acetylsalicylic acid, ASA), dimeric cyclooxygenase 2 (COX2, PTGS2 dimer) can be acetylated. ASA covalently modifies PTGS2 by acetylating a serine residue at position 530 within the cyclooxygenase active site (Lucido et al. 2016). Acetylated PTGS2 dimer (Ac-PTGS2 dimer) changes the oxygenation stereospecificity towards its substrates, perhaps by steric shielding effects (Tosco 2013), producing a shift in lipid mediator production. Ac-PTGS2 dimer expressed in macrophages can be acetylated by ASA, which enables this form to mediate the biosynthesis of precursors of endogenous anti-inflammatory mediators. Ac-PTGS2 dimer is able to incorporate molecular oxygen into \u03c9-3 fatty acid docosahexaenoic acid (DHA), to form 17-hydroxy-docosahexaenoic acid (17-HDHA) (Serhan et al. 2002, Groeger et al. 2010).",
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-9027627.1"},
        "metadata": {
            "source_url": "https://reactome.org/content/detail/R-HSA-9027627.1",
            "source": "Reactome",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "gd_updateTS": "2020-03-05T00:37:18.005Z",
            "creator": "Jassal, B",
            "gd_createTS": "2017-09-05T19:35:03.000Z",
            "creator_orcid": "0000-0002-5039-5405",
        },
    }

    actual = create_nanopub_from_reaction("9027627")
    assert actual

    result = compare(expected, actual)
    assert result


def test_candidateset_subset():

    expected = {
        "assertions": [
            'rxn(reactants(complex(REACTOME:R-HSA-5655280.1!"Activated FGFR4 mutants", p(SP:P19174!PLCG1, pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 1253), pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 472), pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 771), pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 783)), loc(GO:0005886!"plasma membrane"))), products(complex(REACTOME:R-HSA-5655280.1!"Activated FGFR4 mutants", loc(GO:0005886!"plasma membrane")), p(SP:P19174!PLCG1, pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 1253), pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 472), pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 771), pmod(PSIMOD:00048!O4\'-phospho-L-tyrosine, , 783), loc(GO:0005886!"plasma membrane"))))'
        ],
        "annotations": [
            {"type": "Species", "id": "TAX:9606", "label": "Homo sapiens"},
            {"type": "Disease", "id": "DOID:0080006", "label": "bone development disease"},
            {"type": "Disease", "id": "DOID:162", "label": "cancer"},
        ],
        "id": "Reactome_R-HSA-5655336.3",
        "evidence": "Dissociation from the activated receptor quickly follows phosphorylation of PLC-gamma. Phosphorylated PLC-gamma catalyzes the hydrolysis of phosphatidylinositol(4, 5)bisphosphate to generate two second messengers, diacylglycerol and inositol (1,4,5) triphosphate.",
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-5655336.3"},
        "metadata": {
            "source_url": "https://reactome.org/content/detail/R-HSA-5655336.3",
            "source": "Reactome",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "gd_updateTS": "2020-03-31T14:21:45.777Z",
            "creator": "Rothfels, K",
            "gd_createTS": "2012-02-10T01:07:10.000Z",
            "creator_orcid": "0000-0002-0705-7048",
        },
    }

    actual = create_nanopub_from_reaction("5655336")
    assert actual

    result = compare(expected, actual)
    assert result


def test_nestedcomplexes_candidatesets():

    expected = {
        "assertions": [
            'activity(complex(p(REACTOME:R-HSA-8863155.1!"CUL5-box protein"), p(REACTOME:R-HSA-8952553.1!DCUN1Ds), p(REACTOME:R-HSA-8955207.1!COMMDs), p(SP:O60826!CCDC22), p(SP:Q15369!ELOC), p(SP:Q15370!ELOB), p(SP:Q15843!NEDD8, pmod(PSIMOD:00058!N-acetyl-L-methionine, , 1), pmod(PSIMOD:00211!"S-(glycyl)-L-cysteine (Cys-Gly)", , 76)), p(SP:Q93034!CUL5), p(SP:Q969M7!UBE2F, pmod(PSIMOD:00058!N-acetyl-L-methionine, , 1), pmod(PSIMOD:00211!"S-(glycyl)-L-cysteine (Cys-Gly)", , 116)), p(SP:Q9UBF6!RNF7)) directlyIncreases rxn(reactants(complex(p(REACTOME:R-HSA-8863155.1!"CUL5-box protein"), p(REACTOME:R-HSA-8952553.1!DCUN1Ds), p(REACTOME:R-HSA-8955207.1!COMMDs), p(SP:O60826!CCDC22), p(SP:Q15369!ELOC), p(SP:Q15370!ELOB), p(SP:Q15843!NEDD8, pmod(PSIMOD:00058!N-acetyl-L-methionine, , 1), pmod(PSIMOD:00211!"S-(glycyl)-L-cysteine (Cys-Gly)", , 76)), p(SP:Q93034!CUL5), p(SP:Q969M7!UBE2F, pmod(PSIMOD:00058!N-acetyl-L-methionine, , 1), pmod(PSIMOD:00211!"S-(glycyl)-L-cysteine (Cys-Gly)", , 116)), p(SP:Q9UBF6!RNF7)), products(complex(p(REACTOME:R-HSA-8863155.1!"CUL5-box protein"), p(REACTOME:R-HSA-8952553.1!DCUN1Ds), p(REACTOME:R-HSA-8955207.1!COMMDs), p(SP:O60826!CCDC22), p(SP:Q15369!ELOC), p(SP:Q15370!ELOB), p(SP:Q15843!NEDD8, pmod(PSIMOD:00134!N6-glycyl-L-lysine, , 76)), p(SP:Q93034!CUL5, pmod(PSIMOD:01150!"neddylated lysine", , 724)), p(SP:Q9UBF6!RNF7), p(SP:Q969M7!UBE2F, pmod(PSIMOD:00058!N-acetyl-L-methionine, , 1), loc(GO:0005829!cytosol))))'
        ],
        "annotations": [{"type": "Species", "id": "TAX:9606", "label": "Homo sapiens"}],
        "id": "Reactome_R-HSA-8952044.1",
        "evidence": "UBE2F transfers NEDD8 to lysine 724 of CUL5 in the E3 ligase complex (Duda et al, 2008). Neddylation increases the ubiquitination activity of the E3 complex towards its target, and prevents binding of the CUL5 complex with the CAND1 inhibitor (Hori et al, 1999; Duda et al, 2008; Kelsall et al, 2013). Targets of CUL5 RING complexes include a variety of cellular proteins including receptor and non-receptor tyrosine kinases, signaling molecules transcriptional regulators (reviewed in Okamura et al, 2016). CRL5 complexes are also hijacked by viruses such as HIV, HPV and adenovirus among others. Interaction with viral proteins redirects the ubiquitin ligase complex to target host proteins to promote conditions that favor viral propagation (Harada et al, 2002; Mehle et al, 2004; reviewed in Mahon et al, 2014). ",
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-8952044.1"},
        "metadata": {
            "source_url": "https://reactome.org/content/detail/R-HSA-8952044.1",
            "source": "Reactome",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "gd_updateTS": "2020-03-30T20:22:47.942Z",
            "creator": "Rothfels, K",
            "gd_createTS": "2016-12-13T01:49:15.000Z",
            "creator_orcid": "0000-0002-0705-7048",
        },
    }

    actual = create_nanopub_from_reaction("8952044")
    assert actual

    result = compare(expected, actual)
    assert result


def test_nestedcomplexes_entitysets():
    """Result of nested complexes containing entity sets"""
    expected = {
        "annotations": [{"id": "TAX:9606", "label": "Homo sapiens", "type": "Species"}],
        "assertions": [
            'activity(complex(p(REACTOME:R-HSA-4549207.1!"Histone HIST1H2B"), p(REACTOME:R-HSA-68524.3!Ub), p(REACTOME:R-HSA-8963785.1!"Ub-C88-UBE2A,B", pmod(PSIMOD:00211!"S-(glycyl)-L-cysteine (Cys-Gly)")), p(SP:O75150!RNF40), p(SP:Q5VTR2!RNF20), p(SP:Q6P1J9!CDC73), p(SP:Q6PD62!CTR9), p(SP:Q8N7H5!PAF1), p(SP:Q8WVC0!LEO1), p(SP:Q92541!RTF1), p(SP:Q9BTA9!WAC), p(SP:Q9GZS3!WDR61), loc(GO:0005654!nucleoplasm))) directlyIncreases rxn(reactants(complex(p(REACTOME:R-HSA-4549207.1!"Histone HIST1H2B"), p(REACTOME:R-HSA-68524.3!Ub), p(REACTOME:R-HSA-8963785.1!"Ub-C88-UBE2A,B", pmod(PSIMOD:00211!"S-(glycyl)-L-cysteine (Cys-Gly)")), p(SP:O75150!RNF40), p(SP:Q5VTR2!RNF20), p(SP:Q6P1J9!CDC73), p(SP:Q6PD62!CTR9), p(SP:Q8N7H5!PAF1), p(SP:Q8WVC0!LEO1), p(SP:Q92541!RTF1), p(SP:Q9BTA9!WAC), p(SP:Q9GZS3!WDR61), loc(GO:0005654!nucleoplasm))), products(complex(p(SP:O75150!RNF40), p(SP:Q5VTR2!RNF20), loc(GO:0005654!nucleoplasm)), complex(p(SP:Q6P1J9!CDC73), p(SP:Q6PD62!CTR9), p(SP:Q8N7H5!PAF1), p(SP:Q8WVC0!LEO1), p(SP:Q92541!RTF1), p(SP:Q9GZS3!WDR61), loc(GO:0005654!nucleoplasm)), p(REACTOME:R-HSA-6782541.2!"Ub-histone HIST1H2B", pmod(PSIMOD:01148!"ubiquitinylated lysine"), loc(GO:0005654!nucleoplasm)), p(REACTOME:R-HSA-8942138.1!"UBE2A,B", loc(GO:0005654!nucleoplasm)), p(SP:Q9BTA9!WAC, loc(GO:0005654!nucleoplasm))))'
        ],
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-8942101.3"},
        "evidence": "The ubiquitin E3 ligase complex RNF20:RNF40 interacts with the PAF complex (Kim et al. 2009) that is associated with RNA polymerase II via WAC (Zhang and Yu 2011) at transcriptionally active genes (Zhe et al. 2005). RNF20:RNF40 monoubiquitinates nucleosomal histone H2B on lysine-120 (lysine-121 of the unprocessed histone H2B) using UBE2A,B:Ubiquitin as the ubiquitin donor (Zhu et al. 2005, Kim et al. 2009, Zhang and Yu 2011, Zhang et al. 2014, Dickson et al. 2016).  Monoubiquitination of histone H2B leads to methylation of lysine-4 and lysine-79 of histone H3, marks of active chromatin (Zhu et al. 2005). Arsenite binds the RING domains of RNF20 and RNF40 and inhibits the ubiquitination of histone H2B (Zhang et al. 2014).",
        "id": "Reactome_R-HSA-8942101.3",
        "metadata": {
            "creator": "May, B",
            "creator_orcid": "0000-0001-5193-0855",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "source": "Reactome",
            "source_url": "https://reactome.org/content/detail/R-HSA-8942101.3",
        },
    }

    actual = create_nanopub_from_reaction("8942101")
    assert actual

    result = compare(expected, actual)
    assert result


def test_regulator_gene_expression():
    """Regulator of gene expression"""

    expected = {
        "annotations": [{"id": "TAX:9606", "label": "Homo sapiens", "type": "Species"}],
        "assertions": [
            "activity(p(SP:P17861!XBP1, loc(GO:0005654!nucleoplasm))) directlyIncreases p(SP:O14595!CTDSP2, loc(GO:0005654!nucleoplasm))"
        ],
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-1791130.3"},
        "evidence": "The CTDSP2 gene is transcribed to yield mRNA and the mRNA is translated to yield protein.",
        "id": "Reactome_R-HSA-1791130.3",
        "metadata": {
            "creator": "May, B",
            "creator_orcid": "0000-0001-5193-0855",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "source": "Reactome",
            "source_url": "https://reactome.org/content/detail/R-HSA-1791130.3",
        },
    }

    actual = create_nanopub_from_reaction("1791130")
    assert actual

    result = compare(expected, actual)
    assert result


# TODO - fix the MISSING in the object complexes
# def test_protein_int_not_subscriptable():
#     """Protein (int) not subscriptable"""

#     expected = {
#         "annotations": [
#             {
#                 "id": "TAX:9606",
#                 "label": "Homo sapiens",
#                 "type": "Species"
#             }
#         ],
#         "assertions": [
#             "activity(complex(a(CHEBI:29033!\"iron(2+)\"), p(SP:O60568!PLOD3), loc(GO:0005789!\"endoplasmic reticulum membrane\"))) directlyIncreases rxn(reactants(a(CHEBI:18307!UDP-D-galactose, loc(GO:0005788!\"endoplasmic reticulum lumen\")), complex(MISSING(REACTOME:R-HSA-2022990.1!\"Lysyl hydroxylated collagen propeptides\"), a(CHEBI:29033!\"iron(2+)\"), p(SP:O60568!PLOD3)), products(a(CHEBI:17659!UDP, loc(GO:0005788!\"endoplasmic reticulum lumen\")), complex(MISSING(REACTOME:R-HSA-2023001.1!\"Galactosyl-hydroxylysyl collagen propeptides\"), a(CHEBI:29033!\"iron(2+)\"), p(SP:O60568!PLOD3)))"
#         ],
#         "citation": {
#             "uri": "https://reactome.org/content/detail/R-HSA-1981128.2"
#         },
#         "evidence": "The ER membrane-associated enzyme PLOD3 has collagen galactosyltransferase activity (Heikkinen et al. 2000, Wang et al. 2002) though the biological significance of this has been questioned (Schegg et al. 2009).",
#         "id": "Reactome_R-HSA-1981128.2",
#         "metadata": {
#             "creator": "Jupe, S",
#             "creator_orcid": "0000-0001-5807-0069",
#             "license": "CC0",
#             "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
#             "source": "Reactome",
#             "source_url": "https://reactome.org/content/detail/R-HSA-1981128.2"
#         }
#     }

#     actual = create_nanopub_from_reaction("1981128")
#     assert actual

#     result = compare(expected, actual)
#     assert result


def test_other_entity():
    """translocation and other entity example"""

    expected = {
        "annotations": [{"id": "TAX:9606", "label": "Homo sapiens", "type": "Species"}],
        "assertions": [
            'activity(complex(p(REACTOME:R-HSA-1234153.1!HIF-alpha), p(SP:P27540!ARNT), p(SP:Q09472!EP300), p(SP:Q92793!CREBBP), loc(GO:0005654!nucleoplasm))) directlyIncreases p(SP:Q16790!CA9, loc(GO:0005886!"plasma membrane"))'
        ],
        "citation": {"uri": "https://reactome.org/content/detail/R-HSA-1235035.4"},
        "evidence": "The gene encoding carbonic anhydrase IX (CA9) is transcribed to yield mRNA and the mRNA is translated to yield protein. Hypoxia-inducible factor binds the promoter of CA9 and enhances expression of CA9.",
        "id": "Reactome_R-HSA-1235035.4",
        "metadata": {
            "creator": "May, B",
            "creator_orcid": "0000-0001-5193-0855",
            "license": "CC0",
            "license_url": "https://creativecommons.org/publicdomain/zero/1.0",
            "source": "Reactome",
            "source_url": "https://reactome.org/content/detail/R-HSA-1235035.4",
        },
    }

    actual = create_nanopub_from_reaction("1235035")
    assert actual

    result = compare(expected, actual)
    assert result


# Not creating any assertions - but this is a weird Reaction that may not convert to BEL
def test_keyerror_name():
    """Key error - name"""

    expected = {}

    actual = create_nanopub_from_reaction("139952")
    assert actual

    result = compare(expected, actual)
    assert result


def test_recursion_error():
    """Recursion error"""

    expected = {}

    actual = create_nanopub_from_reaction("1235035")
    assert actual

    result = compare(expected, actual)
    assert result


# [info     ] Reaction 72124                 counter=1332 file=./to_bel.py function=convert line=1504 total=12206
# Traceback (most recent call last):
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1133, in get_rna
#     label = doc["referenceEntity"]["geneName"][0]
# KeyError: 'geneName'
# [error    ] 71956	RNA	Error: 'geneName'
#    file=./to_bel.py function=get_rna line=1142
# Traceback (most recent call last):
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1133, in get_rna
#     label = doc["referenceEntity"]["geneName"][0]
# KeyError: 'geneName'
# [error    ] 71911	RNA	Error: 'geneName'


# [info     ] Reaction 9014652               counter=2329 file=./to_bel.py function=convert line=1504 total=12206
# Traceback (most recent call last):
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1052, in get_regulator
#     regulator = process_component(dbid)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1235, in process_component
#     return get_complex(dbid)
#   File "/Users/william/biodati/reactome_to_bel/.venv/lib/python3.7/site-packages/cachetools/func.py", line 74, in wrapper
#     v = func(*args, **kwargs)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 936, in get_complex
#     r = process_component(dbid)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1235, in process_component
#     return get_complex(dbid)
#   File "/Users/william/biodati/reactome_to_bel/.venv/lib/python3.7/site-packages/cachetools/func.py", line 74, in wrapper
#     v = func(*args, **kwargs)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 936, in get_complex
#     r = process_component(dbid)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1235, in process_component
#     return get_complex(dbid)
#   File "/Users/william/biodati/reactome_to_bel/.venv/lib/python3.7/site-packages/cachetools/func.py", line 74, in wrapper
#     v = func(*args, **kwargs)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 936, in get_complex
#     r = process_component(dbid)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1245, in process_component
#     return get_rna(dbid)
#   File "/Users/william/biodati/reactome_to_bel/.venv/lib/python3.7/site-packages/cachetools/func.py", line 74, in wrapper
#     v = func(*args, **kwargs)
#   File "/Users/william/biodati/reactome_to_bel/to_bel.py", line 1144, in get_rna
#     rna = Entity(namespace, id=id_, label=label, stid=stid, stid_version=stid_version, dbid=dbid)
# UnboundLocalError: local variable 'namespace' referenced before assignment
# [error    ] 9014653	Regulator	Error: local variable 'namespace' referenced before assignment
#  file=./to_bel.py function=get_regulator line=1087
# [error    ] Missing assertions             file=./to_bel.py function=convert line=1509 reaction=9014652
