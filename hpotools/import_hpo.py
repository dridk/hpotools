#!/usr/bin/env python
# -*- coding: utf-8 -*-
import networkx as nx
import obonet as obo
import peewee
import pickle
from collections import defaultdict
from operator import itemgetter
import os
from . import model
import hashlib
import datetime
import sys
import logging
from tqdm import tqdm

# =========================================================================================


def read_obo_version(obofile: str) -> str:
    """ read version from obo file """
    with open(obofile, "r") as file:
        for line in file:
            if line.startswith("data-version"):
                return line.split(":")[1].strip()


def dag_from_obo(obofile: str) -> nx.DiGraph:
    """ Return a directed acyclic graph from obo ontology """
    return nx.DiGraph(obo.read_obo(obofile)).reverse()


def tree_from_dag(obograph: nx.DiGraph) -> nx.DiGraph:
    """ Transform a DAG into Tree by splitting node with multiple parents. 
    This function  takes a long time to run  """
    return nx.dag_to_branching(obograph)


def visit_tree(obotree: nx.DiGraph) -> nx.DiGraph:
    """ 
    Visit node tree and assign left and right index according Nested set algorithm 
    https://en.wikipedia.org/wiki/Nested_set_model
    """
    index = 0
    depth = 0
    parent = None
    for i in nx.dfs_labeled_edges(obotree):
        node_name_1 = i[0]
        node_name_2 = i[1]
        sens = i[2]

        if sens == "forward":
            depth += 1
            parent = node_name_1
            obotree.nodes[node_name_2].update(
                {"left": index, "depth": depth, "parent": parent}
            )

        if sens == "reverse":
            depth -= 1
            obotree.nodes[node_name_2].update({"right": index, "depth": depth})

        index += 1

    return obotree


# =========================================================================================
def import_obo_file(obo_filename):

    print("loading obo file ", obo_filename)
    obograph = dag_from_obo(obo_filename)

    print("converting dag to tree. It takes a long time to run ...")
    obotree = tree_from_dag(obograph)

    print("Visiting tree and assign left and right index ")
    tree = visit_tree(obotree)

    print("Create database")

    model.Informations(
        hpo_version=read_obo_version(obo_filename),
        hpo_md5=hashlib.md5(open(obo_filename, "rb").read()).hexdigest(),
        saved_by="Hpo tools",
        version=str(datetime.datetime.now()),
    ).save()

    #  Store all terms
    print("Storing terms")
    all_terms = dict()

    with model.db.transaction() as txn:
        for node in obograph.nodes(data=True):
            term = model.Terms()
            term.hpo = node[0]
            term.name = node[1].get("name")
            term.definition = node[1].get("def", "")
            term.comment = node[1].get("comment", "")

            term.save()

            all_terms[term.hpo] = term

    # Compute SQL
    print("Storing nodes")
    with model.db.transaction() as txn:
        for i in obotree.nodes(data=True):
            item = model.Nodes()
            item.name = i[0]
            item.left = i[1]["left"]
            item.right = i[1]["right"]
            item.depth = i[1]["depth"]
            item.term = all_terms[i[1]["source"]]

            try:
                item.parent = model.Nodes.get(model.Nodes.name == i[1]["parent"])
            except:
                print("cannot store with item.parent=", i[1]["parent"])

            item.save()


# # =========================================================================================


def import_gene_file(gene_filename):
    """ Import gene file 
    http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt 
    """

    def chunks(l, n):
        # For item i in a range that is a length of l,
        for i in range(0, len(l), n):
            # Create an index range for l of n items:
            yield l[i : i + n]

    # save all hpo in memory map
    hpo = {}
    for term in model.Terms.select():
        hpo[term.hpo] = term

    # Load all genes
    genes = set()
    print("create genes tables")
    with open(gene_filename, "r") as file:
        next(file)
        for line in file:
            row = line.rstrip().split("\t")
            # entrezid / gene symbole
            genes.add((int(row[2]), row[3]))

    data_source = [{"entrez_id": i[0], "name": i[1]} for i in genes]

    with model.db.atomic():
        model.Genes.insert_many(data_source).execute()

        # for data_dict in data_source:
        #    Genes.create(**data_dict)

    # save all hpo in memory map
    genes = {}
    for gene in model.Genes.select():
        genes[gene.name] = gene

    # Save hpo_has_gene
    hpo_genes = []
    print("create hpo_has_gene tables")

    with open(gene_filename, "r") as file:
        next(file)
        for line in file:
            row = line.rstrip().split("\t")
            gene = genes[row[3]]
            term = hpo[row[0]]

            hpo_genes.append({"term": term.id, "gene": gene.id})

    with model.db.atomic():
        for items in tqdm(list(chunks(hpo_genes, 100))):
            model.Genes_has_Terms.insert_many(items).execute()


def import_phenotype_file(phenotype_file):
    """ Import phenotype file 
    http://compbio.charite.de/jenkins/job/hpo.annotations.current/lastSuccessfulBuild/artifact/current/phenotype.hpoa
    """

    def chunks(l, n):
        # For item i in a range that is a length of l,
        for i in range(0, len(l), n):
            # Create an index range for l of n items:
            yield l[i : i + n]

    #  disease
    disease_file = phenotype_file

    ids = set()
    diseases = []

    with open(disease_file, "r") as file:
        with model.db.transaction() as txn:
            for line in file:
                if line.startswith("#"):
                    continue
                row = line.rstrip().split("\t")

                database_id = row[0]
                name = row[1]

                if database_id not in ids:
                    ids.add(database_id)
                    diseases.append({"database_id": database_id, "name": name})

    with model.db.atomic():
        model.Diseases.insert_many(diseases).execute()

    # save all hpo in memory map
    hpo = {}
    for term in model.Terms.select():
        hpo[term.hpo] = term

    print(hpo)

    # save all disease in memory map
    diseases = {}
    for disease in model.Diseases.select():
        diseases[disease.database_id] = disease

    print("save diseases relation")
    relations = []
    with open(disease_file, "r") as file:
        with model.db.transaction() as txn:
            for line in file:
                if line.startswith("#"):
                    continue
                row = line.rstrip().split("\t")
                print(row)
                disease_id = diseases[row[0]].id
                hpo_id = hpo[row[3]].id
                qualifier = not (row[2] == "NOT")
                evidence = row[5]
                aspect = row[10]

                relations.append(
                    {
                        "disease": disease_id,
                        "term": hpo_id,
                        "qualifier": qualifier,
                        "evidence": evidence,
                        "aspect": aspect,
                    }
                )

    with model.db.atomic():
        for relation in tqdm(list(chunks(relations, 100))):
            model.Disease_has_Terms.insert_many(relation).execute()

    # Compute diseases frequency
    print("Compute diseases count per terms")

    terms = []
    for term in tqdm(model.Terms.select()):
        count = model.db.execute_sql(
            f"""SELECT  COUNT(DISTINCT d.disease_id) as count FROM nodes
            INNER JOIN (SELECT left, right, term_id FROM nodes WHERE term_id = {term.id}) as root ON nodes.left >= root.left AND nodes.right <= root.right
            INNER JOIN disease_has_terms d ON d.term_id = nodes.term_id"""
        ).fetchone()[0]
        term.disease_count = count
        terms.append(term)

    model.Terms.bulk_update(terms, fields=[model.Terms.disease_count])

    print("done")


def import_all(db_file, obo_file, gene_file, pheno_file):

    model.db.init(db_file)
    model.create_database_shema()

    import_obo_file(obo_file)
    import_gene_file(gene_file)

    import_phenotype_file(pheno_file)
