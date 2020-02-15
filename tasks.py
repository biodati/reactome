from invoke import task

import settings


@task()
def start_arango(c):
    """Start arangodb docker container"""

    c.run(
        "docker run -d --name reactome_arangodb -e ARANGO_NO_AUTH=1 -e ARANGO_STORAGE_ENGINE=rocksdb -p 18529:8529 -v /Users/william/biodati/reactome_to_bel/data/arangodb:/var/lib/arangodb3 arangodb"
    )


@task()
def restart_arango(c):
    """Start arangodb docker container"""

    c.run("docker stop reactome_arangodb; docker rm reactome_arangodb")
    c.run(
        "docker run -d --name reactome_arangodb -e ARANGO_NO_AUTH=1 -e ARANGO_STORAGE_ENGINE=rocksdb -p 18529:8529 -v /Users/william/biodati/reactome_to_bel/data/arangodb:/var/lib/arangodb3 arangodb"
    )


@task()
def stop_arango(c):
    """Stop arangodb docker container"""

    c.run("docker stop reactome_arangodb")
