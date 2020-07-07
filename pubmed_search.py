import argparse
import sqlite3
import os
from Bio import Entrez
from pubmed_lookup import PubMedLookup, Publication
import logging

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
    
    parser = argparse.ArgumentParser(description='Get publication data from PubMed.')
    parser.add_argument(
        'email', type=str, help='Identify yourself so NCBI can contact you in case of excessive requests')
    parser.add_argument('search_query', type=str,
                        help='Search term for Pubmed database')
    parser.add_argument('-s', '--start', type=int,
                        help='Sequential index of the first PMID in the retrieved set. (default: 0)', default=0)
    parser.add_argument('-m', '--max', type=int,
                        help='Total number of PMIDs from the retrieved set. Max: 100,000 records. (default: 100000)', default=100000)
    parser.add_argument('-o', '--output', type=str,
                        help='Specify an existing database. (default: create or connect to a database named after the search query).')                    
    parser.add_argument('-r', '--retry', action='store_true',
                        help='Try to fetch again any PMIDs marked as failed in the database.', default=False)
    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    if args.output:
        database = Database(args.output)
    else:
        database = Database(f'{args.search_query}.db')
    
    search = PubmedSearch(args.email)
    logging.info(f'Getting PubMed IDs (PMID) for articles related to search query "{args.search_query}"...')
    pubmed_ids = search.get_ids(args.search_query, args.start, args.max)
    logging.info(f'{len(pubmed_ids)} PMIDs where found.')
    logging.info('Saving PMIDs to the database...')
    for id in pubmed_ids:
        database.insert_id(id)
    logging.info('PMIDs saved.')