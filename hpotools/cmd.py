import sqlite3


def get_terms(conn: sqlite3.Connection, parent: str, maxdepth=0, showtree=True):

    q = conn.execute(
        f"""SELECT terms.hpo, terms.name, nodes.depth FROM nodes
INNER JOIN (SELECT left, right FROM nodes WHERE term_id = (SELECT id FROM terms WHERE hpo = "{parent}")) as root ON nodes.left >= root.left AND nodes.right <= root.right     
INNER JOIN terms ON terms.id = nodes.term_id"""
    )

    first_depth = None
    for record in q:
        hpo, name, depth = record

        if first_depth is None:
            first_depth = depth

        relative_depth = depth - first_depth

        if relative_depth <= maxdepth or maxdepth is 0:

            if showtree:
                indent = "\t" * relative_depth
                print(indent, "├──", hpo, name)
            else:
                print(hpo, name, sep="\t")


def get_diseases(conn: sqlite3.Connection, terms: list):
    pass
