import argparse
import sqlite3
import os
from Bio import Entrez
from pubmed_lookup import PubMedLookup, Publication

class Database(object):
    def __init__(self, database_name: str):
        self._connection = self.__create_connection(database_name)
        self._create_tables()

    def __create_connection(self, database_name: str) -> sqlite3.Connection:
        database = f'{database_name}.db'
        conn = sqlite3.connect(database)
        conn.row_factory = lambda cursor, row: row[0]
        return conn

    def _create_tables(self):
        try:
            cursor = self._connection.cursor()
            cursor.execute("""PRAGMA foreign_keys=on;""")
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS pmids (
                pmid INTEGER NOT NULL PRIMARY KEY,
                publication_fetched BOOLEAN NOT NULL CHECK (publication_fetched IN (0,1))
            );""")

            cursor.execute("""
            CREATE TABLE IF NOT EXISTS publications (
                pmid INTEGER NOT NULL,
                title TEXT,
                authors TEXT,
                journal TEXT,
                year INTEGER,
                month INTEGER,
                day INTEGER,
                abstract TEXT,
                FOREIGN KEY(pmid) REFERENCES pmids(pmid)
            );""")

            self._connection.commit()
        finally:
            cursor.close()

class PubmedSearch(object):
    def __init__(self, email: str):
        self._email = email

    def get_ids(self, search_query: str, start: int = 0, max: int = 100000) -> list:
        Entrez.email = self._email
        handle = Entrez.esearch(
            db="pubmed", term=search_query, retstart=start, retmax=max)
        record = Entrez.read(handle)
        return list(record["IdList"])

    def get_publication(self, pubmedID: str) -> Publication:
        try:
            lookup = PubMedLookup(pubmedID, self._email)
            return Publication(lookup)
        except TypeError:
            return

if __name__ == 'main':
    pass