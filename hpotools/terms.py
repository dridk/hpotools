import sqlite3


def get_terms(conn: sqlite3.Connection, parent: str, depth=4, show_tree=True):

    q = conn.execute(
        f"""SELECT terms.hpo, terms.name FROM nodes
INNER JOIN (SELECT left, right FROM nodes WHERE term_id = (SELECT id FROM terms WHERE hpo = "{parent}")) as root ON nodes.left > root.left AND nodes.right < root.right     
INNER JOIN terms ON  terms.id = nodes.term_id"""
    )

    for record in q:
        print(record)
