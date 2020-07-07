import argparse
import sqlite3
import os
from Bio import Entrez
from pubmed_lookup import PubMedLookup, Publication

class Database(object):
    def __init__(self, database: str):
        self._connection = self.__create_connection(database)
        self._create_tables()

    def __create_connection(self, database: str) -> sqlite3.Connection:
        conn = sqlite3.connect(database)
        return conn

    def _create_tables(self):
        try:
            cursor = self._connection.cursor()
            cursor.execute("""PRAGMA foreign_keys=on;""")
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS pmids (
                pmid INTEGER NOT NULL PRIMARY KEY,
                new BOOLEAN NOT NULL CHECK (new IN (0,1)),
                failed BOOLEAN NOT NULL CHECK (failed IN (0,1))
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
    
    def insert_id(self, id: int):
        try:
            cursor = self._connection.cursor()
            cursor.execute("INSERT INTO pmids (pmid, new, failed) VALUES (?, ?, ?)", [id, 1, 0])
            self._connection.commit()
        finally:
            cursor.close()
    
    def insert_publication(self, id: int, publication: Publication):
        try:
            cursor = self._connection.cursor()
            cursor.execute("""INSERT INTO publications (pmid, title, authors, journal, year, month, day, abstract)
                                    VALUES (?,?,?,?,?,?,?,?)""",
                                    [id,publication.title, publication.authors, publication.journal,
                                    publication.year, publication.month, publication.day, publication.abstract])
            self._connection.commit()
        finally:
            cursor.close()

    def update_fetched_publication(self, id: int, new:int=0, failed:int=0):
        try:
            cursor = self._connection.cursor()
            cursor.execute("""UPDATE pmids SET new = ?, failed = ? WHERE pmid = ?""", [new, failed, id])
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

if __name__ == '__main__':
    pass