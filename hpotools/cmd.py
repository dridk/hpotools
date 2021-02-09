import sqlite3
import math
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
import vcf


def get_terms(conn: sqlite3.Connection, parent: str, maxdepth=0):
    """ Return child terms from parent terms using a maximum depth """

    yield from conn.execute(
        f"""SELECT DISTINCT terms.hpo, terms.name, nodes.depth FROM nodes
INNER JOIN (SELECT left, right FROM nodes WHERE term_id = (SELECT id FROM terms WHERE hpo = "{parent}")) as root ON nodes.left >= root.left AND nodes.right <= root.right     
INNER JOIN terms ON terms.id = nodes.term_id WHERE nodes.depth <= {maxdepth}"""
    )

    # first_depth = None
    # for record in q:
    #     hpo, name, depth = record

    #     if first_depth is None:
    #         first_depth = depth

    #     relative_depth = depth - first_depth

    #     if relative_depth <= maxdepth or maxdepth is 0:

    #         if showtree:
    #             indent = "\t" * relative_depth
    #             print(indent, "├──", hpo, name)
    #         else:
    #             print(hpo, name, sep="\t")


def get_terms_by_name(conn: sqlite3.Connection, query: str):
    """ Select hpo term by a name. The query use a SQL LIKE expression  """
    yield from conn.execute(
        f"""SELECT DISTINCT terms.hpo, terms.name FROM terms WHERE terms.name LIKE '%{query}%' """
    )


def get_term_by_nodeid(conn: sqlite3.Connection, node_id: int):
    return conn.execute(
        f"""SELECT terms.hpo, terms.name FROM terms, nodes WHERE nodes.term_id = terms.id and nodes.id = {node_id} """
    ).fetchone()


def get_term_by_id(conn: sqlite3.Connection, term_id: int):
    return conn.execute(
        f"SELECT terms.hpo, terms.name FROM terms WHERE terms.id = {term_id}"
    ).fetchone()


def get_genes(conn: sqlite3.Connection, term: str):
    for record in conn.execute(
        f"""
        SELECT DISTINCT genes.name FROM terms , genes
        INNER JOIN genes_has_terms ON genes_has_terms.gene_id = genes.id AND genes_has_terms.term_id = terms.id
        WHERE terms.hpo = "{term}"
        ORDER BY genes.name """
    ):
        yield record[0]


def get_genes_by_terms(conn: sqlite3.Connection, terms: list):
    """ Return genes from terms"""

    setlist = []
    for term in terms:
        setlist += list(get_genes(conn, term))

    return setlist


def get_diseases(conn: sqlite3.Connection, term: str):
    yield from conn.execute(
        f"""SELECT diseases.database_id,  diseases.name FROM diseases
        INNER JOIN terms ON disease_has_terms.term_id = terms.id 
        INNER JOIN disease_has_terms ON  disease_has_terms.disease_id = diseases.id
        WHERE terms.hpo = "{term}"
        """
    )


def get_diseases_by_terms(conn: sqlite3.Connection, terms: list):
    """ Return genes from terms"""

    setlist = []
    for term in terms:
        setlist += list(get_diseases(conn, term))

    return setlist


def get_common_ancestors(conn: sqlite3.Connection, term1: str, term2: str):
    """Get term of the common ancestor between two terms
    
    Args:
        conn (sqlite3.Connection): sqlite3 Connection
        term1 (str): HPO Term1
        term2 (str): HPO term2 
    """

    q = "SELECT nodes.term_id, MIN(left) as min_left, MIN(right) as min_right FROM nodes WHERE nodes.term_id IN (SELECT terms.id FROM terms WHERE hpo IN ('{}', '{}')) GROUP BY term_id".format(
        term1, term2
    )

    # TODO... Actually I returned the first common ancestor between duplicated nodes
    term_id = []
    for record in conn.execute(q):
        _, min_left, max_right = record
        q = f"SELECT nodes.term_id FROM nodes WHERE nodes.left < {min_left}  AND nodes.right > {max_right} ORDER BY left DESC LIMIT 1 "
        (term_id,) = conn.execute(q).fetchone()
        return term_id


def get_term_entropy(conn: sqlite3.Connection, term: str):
    """Compute the entropy score of a terms 
    
    Args:
        conn (sqlite3 / Connection): sqlite3.Connection
        term (str): HPO Terms 
    """

    # Get total diseases count   => correspond to disease count of root node
    (total,) = conn.execute(
        "SELECT disease_count FROM terms WHERE hpo = 'HP:0000001'"
    ).fetchone()

    # Get disease count for the specific term
    q = f"SELECT disease_count FROM terms where hpo = '{term}'"
    (result,) = conn.execute(q).fetchone()

    p = float(result) / float(total)

    entropy = -math.log2(p)
    return entropy


def get_resnik_score(conn: sqlite3.Connection, term1: str, term2: str):
    """Compute resnik score similarity between 2 hpo terms 
    
    Args:
        conn (sqlite3.Connection): sqlite3.Connection
        term1 (str): hpo term
        term2 (str): hpo term
    """

    term_id = get_common_ancestors(conn, term1, term2)
    (lca,) = conn.execute(f"SELECT hpo FROM terms WHERE id = {term_id}").fetchone()

    return get_term_entropy(conn, lca)


def get_term_similarity(
    conn: sqlite3.Connection, term1: str, term2: str, metrics="resnik"
):
    if metrics == "resnik":
        return get_resnik_score(conn, term1, term2)


def get_phenotype_similarity(
    conn: sqlite3.Connection, hpo_set1: set(), hpo_set2: set(), metrics="resnik"
):

    results = []
    for i in hpo_set1:
        for j in hpo_set2:
            results.append(get_term_similarity(conn, i, j, metrics))

    return sum(results) / len(results)


def get_pairwise_matrix(conn: sqlite3.Connection, pheno: dict, metrics="resnik"):

    phenotypes = list(pheno.keys())
    X = np.array(phenotypes).reshape(-1, 1)

    def similarity(key1, key2):

        set1 = pheno[key1[0]]
        set2 = pheno[key2[0]]
        return get_phenotype_similarity(conn, set1, set2, metrics)

    dist = pdist(X, metric=similarity)
    df = pd.DataFrame(squareform(dist))
    df.columns = phenotypes
    df.index = phenotypes

    return df


def annotate_vcf(
    conn: sqlite3.Connection,
    input_file: str,
    output_file: str,
    max_freq=0.01,
    use_name=False,
):

    # get total diseases
    total_disease = int(conn.execute("SELECT max(id) FROM diseases").fetchone()[0])

    reader = vcf.Reader(filename=input_file)
    writer = vcf.Writer(open(output_file, "w"), reader)

    for record in reader:
        genes = set()
        if "ANN" in record.INFO:
            for sub_record in record.INFO["ANN"]:
                fields = sub_record.split("|")
                # gene field is in 3 column :
                genes.add(fields[3])

        # query hpo terms related to genes

        sql_field = "terms.name" if use_name else "terms.hpo"
        gene_filter = "(" + ",".join(f"'{g}'" for g in genes) + ")"

        sql = f"""SELECT {sql_field}, CAST(terms.disease_count AS FLOAT) / {total_disease} as freq  FROM genes
LEFT JOIN genes_has_terms ON genes.id = genes_has_terms.gene_id
LEFT JOIN terms ON genes_has_terms.term_id = terms.id
WHERE genes.name IN {gene_filter} AND freq < {max_freq}"""

        terms = []
        for sql_record in conn.execute(sql):
            terms.append(sql_record[0])

        # write hpo terms
        record.INFO["HpoTools"] = ";".join(terms)

        writer.write_record(record)


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3864022/
