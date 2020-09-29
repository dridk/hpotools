import sqlite3
import math
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd 

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


def get_diseases(conn: sqlite3.Connection, term: str):
    q = conn.execute(
        f"""SELECT diseases.database_id,  diseases.name FROM terms , diseases
INNER JOIN disease_has_terms ON disease_has_terms.term_id = terms.id and disease_has_terms.disease_id = diseases.id

WHERE terms.hpo = "{term}"
        """
    )

    for record in q:
        database, name = record
        print(database, name, sep="\t")



def get_common_ancestors(conn: sqlite3.Connection, term1:str, term2: str):
    """Get term of the common ancestor between two terms
    
    Args:
        conn (sqlite3.Connection): sqlite3 Connection
        term1 (str): HPO Term1
        term2 (str): HPO term2 
    """

    q = "SELECT nodes.term_id, MIN(left) as min_left, MIN(right) as min_right FROM nodes WHERE nodes.term_id IN (SELECT terms.id FROM terms WHERE hpo IN ('{}', '{}')) GROUP BY term_id".format(term1, term2)

    # TODO... Actually I returned the first common ancestor between duplicated nodes 
    term_id = []
    for record in conn.execute(q):
        _, min_left, max_right = record 
        q = f"SELECT nodes.term_id FROM nodes WHERE nodes.left < {min_left}  AND nodes.right > {max_right} ORDER BY left DESC LIMIT 1 "
        term_id, = conn.execute(q).fetchone()
        return term_id


def get_term_entropy(conn: sqlite3.Connection, term:str):
    """Compute the entropy score of a terms 
    
    Args:
        conn (sqlite3 / Connection): sqlite3.Connection
        term (str): HPO Terms 
    """

    # Get total diseases count   => correspond to disease count of root node 
    total, = conn.execute("SELECT disease_count FROM terms WHERE hpo = 'HP:0000001'").fetchone()
    
    # Get disease count for the specific term 
    q = f"SELECT disease_count FROM terms where hpo = '{term}'"
    result, = conn.execute(q).fetchone()

    p = float(result) / float(total)
    
    entropy = -math.log2(p)
    return entropy

def get_resnik_score(conn: sqlite3.Connection, term1:str, term2:str):
    """Compute resnik score similarity between 2 hpo terms 
    
    Args:
        conn (sqlite3.Connection): sqlite3.Connection
        term1 (str): hpo term
        term2 (str): hpo term
    """
    
    term_id = get_common_ancestors(conn, term1, term2)
    lca, = conn.execute(f"SELECT hpo FROM terms WHERE id = {term_id}").fetchone()

    return get_term_entropy(conn, lca)



def get_term_similarity(conn: sqlite3.Connection, term1:str, term2:str, metrics = "resnik"):
    if metrics == "resnik":
        return get_resnik_score(conn, term1, term2)

def get_phenotype_similarity(conn: sqlite3.Connection, hpo_set1 : set(), hpo_set2: set(), metrics= "resnik"):
    
    results = []
    for i in hpo_set1:
        for j in hpo_set2:
            results.append(get_term_similarity(conn, i, j, metrics))

    return sum(results)/len(results)


def get_pairwise_matrix(conn: sqlite3.Connection, pheno: dict , metrics="resnik"):

    phenotypes = list(pheno.keys())
    X = np.array(phenotypes).reshape(-1,1)

    def similarity(key1, key2):
        
        set1 = pheno[key1[0]]
        set2 = pheno[key2[0]]
        return get_phenotype_similarity(conn, set1, set2, metrics)

    dist = pdist(X, metric = similarity)
    df = pd.DataFrame(squareform(dist))
    df.columns = phenotypes 
    df.index = phenotypes

    return df  

   






# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3864022/



