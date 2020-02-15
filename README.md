# Reactome conversion notes

## Complexes and EntitySets

-   Complexes are flattened and only non-redundant members are shown since in BEL we don't track stochiometry
-   EntitySets are treated like Protein families - both DefinedSets and CandidateSets are treated the same.

## Citations

-   The URL for the Reactome Reaction is used for the citation. There may be multiple literature references supporting a Reaction.

## Metadata

-   gd_updateTS is the date that the nanopub was created
-   gd_createTS is the data the Reactome reaction record was created or updated

## Arangodb notes

data volume in docker: /var/lib/arangodb3

    docker run -e ARANGO\_NO\_AUTH=1 -e ARANGO\_STORAGE\_ENGINE=rocksdb -p 18529:8529 -v /Users/william/biodati/reactome\_to\_bel/data/arangodb:/var/lib/arangodb3 arangodb
