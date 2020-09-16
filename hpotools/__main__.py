import argparse
import sys
import os
import urllib.request
import sqlite3
from .cmd import *

HPO_URL = "https://github.com/dridk/hpo2sqlite/releases/download/latest/hpo.db"


class HpoParser(object):
    def __init__(self):

        parser = argparse.ArgumentParser(
            description="Hpo tools",
            usage="""hpotools <command> [<args>] 
subcommand are : 

    init   Init sqlite database path
    terms  query hpo terms 
    search search a terms according his name
    gene   get gene of a terms 
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

    def _get_conn(self):
        dbfile = os.path.join(os.path.expanduser("~"), ".hpotools", "hpo.db")
        if not os.path.exists(dbfile):
            print("Error: hpo.db doesn't exists. Please use hpotools init")
            exit(1)

        return sqlite3.connect(dbfile)


if __name__ == "__main__":

    HpoParser()
