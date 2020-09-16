import sqlite3


def get_terms(conn: sqlite3.Connection, parent: str, maxdepth=0, showtree=True):

    q = conn.execute(
        f"""SELECT DISTINCT terms.hpo, terms.name, nodes.depth FROM nodes
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


def get_terms_by_desc(conn: sqlite3.Connection, query: str):
    q = conn.execute(
        f"""SELECT DISTINCT terms.hpo, terms.name FROM terms WHERE terms.name LIKE '%{query}%' """
    )

    for record in q:
        hpo, name = record
        print(hpo, name, sep="\t")


def get_genes(conn: sqlite3.Connection, term: str):
    q = conn.execute(
        f"""
SELECT DISTINCT genes.name FROM terms , genes
INNER JOIN genes_has_terms ON genes_has_terms.gene_id = genes.id AND genes_has_terms.term_id = terms.id
WHERE terms.hpo = "{term}"
ORDER BY genes.name

        """
    )

    for record in q:
        print(record[0])


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3864022/
