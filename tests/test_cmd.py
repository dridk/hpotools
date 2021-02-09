import pytest

from hpotools import cmd
from hpotools.settings import HPOTOOLS_DIR, HPO_URL
import urllib.request
import os
import sqlite3
import pandas as pd

ROOT_TERM = "HP:0000001"
ROOT_CHILDREN = [
    ("HP:0000118", "Phenotypic abnormality"),
    ("HP:0000005", "Mode of inheritance"),
    ("HP:0031797", "Clinical course"),
    ("HP:0012823", "Clinical modifier"),
    ("HP:0032443", "Past medical history"),
    ("HP:0032223", "Blood group"),
    ("HP:0040279", "Frequency"),
]


@pytest.fixture(scope="session", autouse=True)
def conn():

    print("Download")
    if not os.path.exists(HPOTOOLS_DIR):
        os.mkdir(HPOTOOLS_DIR)

    dbfile = os.path.join(HPOTOOLS_DIR, "hpo.db")
    print(dbfile)
    if not os.path.exists(dbfile):
        urllib.request.urlretrieve(HPO_URL, dbfile)

    return sqlite3.connect(dbfile)


def test_get_terms(conn):

    # Get all terms from first level with depth = 1
    for term, name, level in cmd.get_terms(conn, parent=ROOT_TERM, maxdepth=1):
        assert level <= 1
        if term != ROOT_TERM:
            assert term in [i[0] for i in ROOT_CHILDREN]


def test_get_terms_by_nodeid(conn):

    assert cmd.get_term_by_nodeid(conn, 1) == ("HP:0000001", "All")


def test_get_terms_by_id(conn):
    assert cmd.get_term_by_id(conn, 1) == ("HP:0000001", "All")


def test_get_terms_by_name(conn):

    word = "headache"
    for hpo, name in cmd.get_terms_by_name(conn, word):
        assert word in name.lower()


def test_get_genes(conn):
    # Look gene of sclerodactylye
    excepted = ["LBR", "RSPO1", "SMARCAD1"]
    for gene in cmd.get_genes(conn, "HP:0011838"):
        print(gene)
        assert gene in excepted


def test_get_genes_by_terms(conn):

    # There is only one gene RSPO1 in common between HP:0007410 and HP:0011838
    common_gene = cmd.get_genes_by_terms(conn, ("HP:0007410", "HP:0011838"))

    serie = pd.Series(common_gene).value_counts()

    assert len(serie[serie == 2])
    assert serie[serie == 2].index[0] == "RSPO1"


def test_get_diseases(conn):

    assert "OMIM:257980" in [
        record[0] for record in cmd.get_diseases(conn, "HP:0001810")
    ]


def test_get_diseases_by_terms(conn):
    diseases = cmd.get_diseases_by_terms(conn, ("HP:0007410", "HP:0011838"))


def test_common_ancestor(conn):

    assert cmd.get_common_ancestors(conn, ROOT_CHILDREN[0][0], ROOT_CHILDREN[2][0]) == 1

    nodeid = cmd.get_common_ancestors(conn, "HP:0000234", "HP:0000464")
    assert cmd.get_term_by_id(conn, nodeid)[0] == "HP:0000152"


def test_get_entropy(conn):
    cmd.get_term_entropy(conn, "HP:0000001") == 0

    root_term = cmd.get_term_entropy(conn, "HP:0000769")
    leaf_term = cmd.get_term_entropy(conn, "HP:0100013")

    assert root_term < leaf_term
