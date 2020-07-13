import argparse
import sqlite3
import os
from Bio import Entrez
from pubmed_lookup import PubMedLookup, Publication
import logging
import sys

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
    
    def get_new_ids(self):
        try:
            cursor = self._connection.cursor()
            cursor.execute("""SELECT pmid FROM pmids WHERE new = 1""")
            self._connection.commit()
            result = cursor.fetchall()
            cursor.close()
            ids = [id[0] for id in result]
            return ids
        finally:
            cursor.close()
    
    def get_failed_ids(self):
        try:
            cursor = self._connection.cursor()
            cursor.execute("""SELECT pmid FROM pmids WHERE failed = 1""")
            self._connection.commit()
            result = cursor.fetchall()
            cursor.close()
            ids = [id[0] for id in result]
            return ids
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
        lookup = PubMedLookup(pubmedID, self._email)
        return Publication(lookup)

def _parse_args(args: list):
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
    return parser.parse_args()

if __name__ == '__main__':
    
    args = _parse_args(sys.argv[1:])

    logging.basicConfig(filename=f'{args.search_query}.log', format='%(asctime)s - %(message)s', level=logging.INFO)
    logging.info(f'Search query: {args.search_query}')

    if args.output:
        database = Database(args.output)
    else:
        database = Database(f'{args.search_query}.db')

    print(f'Getting PubMed IDs (PMID) for articles related to search query "{args.search_query}"...')
    search = PubmedSearch(args.email)
    pubmed_ids = search.get_ids(args.search_query, args.start, args.max)
    print(f'{len(pubmed_ids)} PMIDs where found.')
    logging.info(f'{len(pubmed_ids)} PMIDs where found.')
    
    print('Saving PMIDs to the database...')
    for id in pubmed_ids:
        try:
            database.insert_id(id)
        except sqlite3.IntegrityError:
            continue

    if args.retry:
        pubmed_ids = database.get_new_ids() + database.get_failed_ids()
        print(f'{len(database.get_new_ids())} PMIDs are marked as NEW in the database.')
        print(f'{len(database.get_failed_ids())} PMIDs are marked as FAILED in the database.')
    else:
        pubmed_ids = database.get_new_ids()
        print(f'{len(pubmed_ids)} PMIDs are marked as NEW in the database.')

    print('Fetching publications...\n')
    count_success = 0
    count_failed = 0
    total = len(pubmed_ids)

    for id in pubmed_ids:
        try:
            publication = search.get_publication(id)
            database.insert_publication(id, publication)
            database.update_fetched_publication(id)
            count_success = count_success + 1
        except (TimeoutError, RuntimeError, TypeError):
            logging.exception(f'ID: {id}')
            database.update_fetched_publication(id, failed=1)
            count_failed = count_failed + 1
        finally:
            print(f'  Successful: {count_success} | Failed: {count_failed} | Remaining: {total-count_success-count_failed}', end='\r', flush=True)
    
    print("\nFinished!")
    logging.info(f'Successful: {count_success} | Failed: {count_failed}')