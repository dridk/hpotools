import networkx as nx
from hpotools import import_hpo as ii
from hpotools import model
import tempfile
import sqlite3

OBOFILE = "examples/hptest.obo"
DBFILE = tempfile.mkstemp(suffix=".db")[1]
GENEFILE = "examples/gene_test.txt"
PHENOTYPEFILE = "examples/phenotype_test.hpoa"


OBOFILE_TERMS = [
    "HP:0000001",
    "HP:0000005",
    "HP:0000118",
    "HP:0012823",
    "HP:0031797",
    "HP:0032223",
    "HP:0032443",
    "HP:0040279",
]


def test_import_function():

    # test version
    assert ii.read_obo_version(OBOFILE) == "hp/releases/2020-08-11"

    # create graph
    graph = ii.dag_from_obo(OBOFILE)

    # test graph content
    assert list(graph.nodes) == OBOFILE_TERMS

    # Tree from dag
    tree = ii.tree_from_dag(graph)

    assert nx.is_tree(tree)

    final_tree = ii.visit_tree(tree)

    index = 0
    print("\n")
    for node, attr in final_tree.nodes(data=True):

        name = attr["source"]

        if name == OBOFILE_TERMS[0]:

            assert attr["depth"] == 0
            assert attr["right"] - attr["left"] == len(OBOFILE_TERMS) * 2 - 1

        else:
            attr["depth"] = 1
            assert attr["right"] - attr["left"] == 1
        index += 1


def test_import_obo_file():

    model.db.init(DBFILE)
    model.create_database_shema()
    ii.import_obo_file(OBOFILE)

    conn = sqlite3.connect(DBFILE)
    conn.row_factory = sqlite3.Row
    record = conn.execute("SELECT * FROM Informations").fetchone()
    assert record != None

    terms = [record[1] for record in conn.execute("SELECT * FROM terms")]
    assert terms == OBOFILE_TERMS

    for record in conn.execute("SELECT * FROM nodes"):

        record = dict(record)

        if record["id"] == 1:
            assert record["right"] - record["left"] == len(OBOFILE_TERMS) * 2 - 1
            assert record["depth"] == 0
        else:
            assert record["right"] - record["left"] == 1
            assert record["depth"] == 1


def test_import_gene_file():

    model.db.init(DBFILE)
    model.create_database_shema()
    ii.import_obo_file(OBOFILE)

    ii.import_gene_file(GENEFILE)

    conn = sqlite3.connect(DBFILE)
    conn.row_factory = sqlite3.Row

    expected_genes = []
    with open(GENEFILE, "r") as file:
        for line in file:
            line = line.split("\t")
            expected_genes.append(line[3].strip())

    observered_genes = [
        record["name"] for record in conn.execute("SELECT * FROM genes")
    ]

    assert sorted(set(expected_genes)) == sorted(set(observered_genes))


def test_import_phenotype_file():

    model.db.init(DBFILE)
    model.create_database_shema()

    ii.import_obo_file(OBOFILE)
    ii.import_gene_file(GENEFILE)

    ii.import_phenotype_file(PHENOTYPEFILE)

    conn = sqlite3.connect(DBFILE)
    conn.row_factory = sqlite3.Row

    expected_diseases = []
    with open(PHENOTYPEFILE) as file:
        for line in file:
            expected_diseases.append(line.split("\t")[4].strip())

    observed_diseases = [
        record["database_id"] for record in conn.execute("SELECT * FROM diseases")
    ]

    assert expected_diseases == observed_diseases


def test_import_all():

    dbfile = tempfile.mkstemp(suffix=".db")[1]

    ii.import_all(dbfile, OBOFILE, GENEFILE, PHENOTYPEFILE)
