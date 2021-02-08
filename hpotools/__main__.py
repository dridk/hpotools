import argparse
import sys
import os
import urllib.request
import sqlite3
import json
from .cmd import *

HPO_URL = "https://github.com/dridk/hpo2sqlite/releases/download/latest/hpo.db"


class HpoParser(object):
    def __init__(self):

        parser = argparse.ArgumentParser(
            description="Hpo tools",
            usage="""hpotools <command> [<args>]
Examples:

    hpotools init  # Download a precomputed database 
    hpotools search "headache" # Search a terms by name 


subcommand are : 
    init             Init sqlite database path
    build            Build a database from obo file and annotation
    terms            query hpo terms 
    search           search a terms according his name
    gene             get gene of a terms 
    lca              get last common ancestor between two terms 
    term_similarity  Compute similarity between two terms  
    phen_similarity  Compute similarity between two phenotypes
    matrix           Compute similarity matrix between phenotypes
    annotate         Annotate a vcf file previously annotated with snpeff or VEP


            """,
        )

        parser.add_argument("command", help="subcommand to run")

        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command + "_cmd"):
            print("Unrecognized command")
            parser.print_help()
            exit(1)

        getattr(self, args.command + "_cmd")(sys.argv[2:])

    def init_cmd(self, argv):
        parser = argparse.ArgumentParser(description="Init database path")
        parser.add_argument("--database", "-d", type=argparse.FileType("r"))
        args = parser.parse_args(argv)

        dirpath = os.path.join(os.path.expanduser("~"), ".hpotools")
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)

        if args.database is None:
            print("Downloading ", HPO_URL)
            urllib.request.urlretrieve(HPO_URL, os.path.join(dirpath, "hpo.db"))

    def build_cmd(self, argv):
        parser = argparse.ArgumentParser(
            description="""
            Build a sqlite database from hpo obo file and annotations.

            Note: This action can take a while. 
            it is better to download the precomputed database from hpotools init.

    
            """
        )
        parser.add_argument(
            "--hpo",
            "-i",
            type=argparse.FileType("r"),
            help="http://purl.obolibrary.org/obo/hp.obo ",
        )
        parser.add_argument(
            "--phenotype",
            "-p",
            type=argparse.FileType("r"),
            help="http://compbio.charite.de/jenkins/job/hpo.annotations.current/lastSuccessfulBuild/artifact/current/phenotype.hpoa",
        )
        parser.add_argument(
            "--gene",
            "-g",
            type=argparse.FileType("r"),
            help="http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt",
        )
        args = parser.parse_args(argv)

        print(args)

    def terms_cmd(self, argv):
        parser = argparse.ArgumentParser(description="Query terms")
        parser.add_argument("--parent", "-p", type=str)
        parser.add_argument("--showtree", "-t", default=False, action="store_true")
        parser.add_argument("--maxdepth", "-d", type=int, default=0)
        args = parser.parse_args(argv)

        get_terms(
            conn=self._get_conn(),
            parent=args.parent,
            showtree=args.showtree,
            maxdepth=args.maxdepth,
        )

    def search_cmd(self, argv):
        parser = argparse.ArgumentParser(description="Search a  terms")
        parser.add_argument("query", type=str)
        args = parser.parse_args(argv)

        if args.query is not None:
            get_terms_by_desc(self._get_conn(), args.query)

    def gene_cmd(self, argv):
        parser = argparse.ArgumentParser(description="show genes of an HPO term")
        parser.add_argument("term", type=str)
        args = parser.parse_args(argv)

        if args.term is not None:
            get_genes(self._get_conn(), args.term)

    def disease_cmd(self, argv):
        parser = argparse.ArgumentParser(description="show diseases of an HPO term")
        parser.add_argument("term", type=str)
        args = parser.parse_args(argv)

        if args.term is not None:
            get_diseases(self._get_conn(), args.term)

    def lca_cmd(self, argv):
        parser = argparse.ArgumentParser(
            description="get last common ancestors between two hpo terms"
        )
        parser.add_argument("first", type=str)
        parser.add_argument("second", type=str)
        args = parser.parse_args(argv)

        if args.first is not None and args.second is not None:
            conn = self._get_conn()
            term_id = get_common_ancestors(conn, args.first, args.second)
            result = conn.execute(
                f"SELECT terms.hpo, terms.name FROM terms WHERE id = {term_id}"
            ).fetchone()
            print(*result, sep="\t")

    def term_similarity_cmd(self, argv):
        parser = argparse.ArgumentParser(
            description="Compute similarity between two hpo terms "
        )
        parser.add_argument("first", type=str)
        parser.add_argument("second", type=str)
        parser.add_argument(
            "-m", "--metrics", type=str, choices=["resnik"], default="resnik"
        )
        args = parser.parse_args(argv)

        if args.metrics == "resnik":
            print(
                "Resnik score:",
                get_resnik_score(self._get_conn(), args.first, args.second),
            )

    def phen_similarity_cmd(self, argv):
        parser = argparse.ArgumentParser(
            description="Compute similarity between two phenotypes "
        )
        parser.add_argument("-a", "--first", type=str)
        parser.add_argument("-b", "--second", type=str)
        parser.add_argument(
            "-m", "--metrics", type=str, choices=["resnik"], default="resnik"
        )
        args = parser.parse_args(argv)

        print(
            "Score:",
            get_phenotype_similarity(
                self._get_conn(),
                set(args.first.split(",")),
                set(args.second.split(",")),
                metrics=args.metrics,
            ),
        )

        # print("Resnik score:" , get_term_similarity(self._get_conn(), args.first, args.second, metrics = args.metrics))

    def matrix_cmd(self, argv):
        parser = argparse.ArgumentParser(
            description="Compute similarity matrix between two sets of words "
        )
        parser.add_argument("file", type=argparse.FileType("r"))
        parser.add_argument(
            "-m", "--metrics", type=str, choices=["resnik"], default="resnik"
        )
        args = parser.parse_args(argv)

        if args.metrics == "resnik":
            data = json.load(args.file)
            print(get_pairwise_matrix(self._get_conn(), data, args.metrics))
            # print("Resnik score:" , get_resnik_score(self._get_conn(), args.first, args.second))

    def annotate_cmd(self, argv):
        parser = argparse.ArgumentParser(
            description="Annotate a VCF file previously annotated with snpEff or VEP"
        )
        parser.add_argument("-i", "--input", help="annotated vcf file input", type=str)
        parser.add_argument(
            "-o", "--output", type=str, help="annotated vcf file output"
        )
        parser.add_argument(
            "-f",
            "--max-freq",
            type=str,
            help="diseases frequency. A low frequence mean more specificity",
            default=0.01,
        )
        parser.add_argument(
            "-n",
            "--use-name",
            action="store_true",
            help="use HPO name instead of HPO identifier",
            default=False,
        )
        args = parser.parse_args(argv)
        print(args)
        annotate_vcf(
            self._get_conn(),
            args.input,
            args.output,
            max_freq=args.max_freq,
            use_name=args.use_name,
        )

    def _get_conn(self):
        dbfile = os.path.join(os.path.expanduser("~"), ".hpotools", "hpo.db")
        if not os.path.exists(dbfile):
            print("Error: hpo.db doesn't exists. Please use hpotools init")
            exit(1)

        return sqlite3.connect(dbfile)


if __name__ == "__main__":

    HpoParser()
