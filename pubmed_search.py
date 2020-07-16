import argparse
import sqlite3
import os
from Bio import Entrez
from pubmed_lookup import PubMedLookup, Publication
import logging
import sys

class Database():
    def __init__(self, database: str):
        self.__connection = sqlite3.connect(database)
        self.__create_tables()

    @property
    def connection(self) -> sqlite3.Connection:
        return self.__connection

    def __execute_query(self, *args) -> list:
        try:
            cursor = self.connection.cursor()
            cursor.execute(*args)
            self.connection.commit()
            result = cursor.fetchall()
        finally:
            cursor.close()
        return result

    def __create_tables(self):
        self.__execute_query("""PRAGMA foreign_keys=on;""")
        self.__execute_query("""
            CREATE TABLE IF NOT EXISTS pmids (
            pmid INTEGER NOT NULL PRIMARY KEY,
            new BOOLEAN NOT NULL CHECK (new IN (0,1)),
            failed BOOLEAN NOT NULL CHECK (failed IN (0,1))
        );""")

        self.__execute_query("""
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

    def insert_id(self, id: int):
        self.__execute_query("INSERT INTO pmids (pmid, new, failed) VALUES (?, ?, ?)", [id, 1, 0])

    def insert_publication(self, id: int, publication: Publication):
        self.__execute_query("""INSERT INTO publications (pmid, title, authors, journal, year, month, day, abstract)
                                    VALUES (?,?,?,?,?,?,?,?)""",
                             [id, publication.title, publication.authors, publication.journal,
                              publication.year, publication.month, publication.day, publication.abstract])

    def update_fetched_publication(self, id: int, new: int = 0, failed: int = 0):
        self.__execute_query("""UPDATE pmids SET new = ?, failed = ? WHERE pmid = ?""", [new, failed, id])

    def select_new_ids(self) -> list:
        result = self.__execute_query("""SELECT pmid FROM pmids WHERE new = 1""")
        ids = [id[0] for id in result]
        return ids

    def select_failed_ids(self) -> list:
        result = self.__execute_query("""SELECT pmid FROM pmids WHERE failed = 1""")
        ids = [id[0] for id in result]
        return ids

class PubmedSearch():
    def __init__(self, email: str):
        self.__email = email
    
    @property
    def email(self) -> str:
        return self.__email

    @email.setter
    def email(self, email: str):
        self.__email = email

    def get_ids(self, search_query: str, start: int = 0, max: int = 100000) -> list:
        Entrez.email = self.email
        handle = Entrez.esearch(db="pubmed", term=search_query, retstart=start, retmax=max)
        record = Entrez.read(handle)
        return list(record["IdList"])

    def get_publication(self, pubmedID: int) -> Publication:
        lookup = PubMedLookup(pubmedID, self.email)
        return Publication(lookup)

def _parse_args(args: list):
    parser = argparse.ArgumentParser(
        description='Get publication data from PubMed.')
    parser.add_argument(
        'email', type=str, help='Identify yourself so NCBI can contact you in case of excessive requests')
    parser.add_argument('search_query', type=str,
                        help='Search term for Pubmed database')
    parser.add_argument('-s', '--start', type=int,
                        help='Sequential index of the first PMID from retrieved set. (default: 0)', default=0)
    parser.add_argument('-m', '--max', type=int,
                        help='Total number of PMIDs from the retrieved set. Max: 100,000 records. (default: 100000)', default=100000)
    parser.add_argument('-o', '--output', type=str,
                        help='Specify an existing database. (default: create or connect to a database using the search query as a file name).')
    parser.add_argument('-r', '--retry', action='store_true',
                        help='Try to fetch again any PMIDs marked as failed in the database.', default=False)
    return parser.parse_args()


if __name__ == '__main__':
    
    args = _parse_args(sys.argv[1:])
    logging.basicConfig(filename=f'{args.search_query}.log',format='%(asctime)s - %(message)s', level=logging.INFO)
    logging.info(f'Search query: {args.search_query}')

    if args.output:
        database = Database(args.output)
    else:
        database = Database(f'{args.search_query}.db')

    print(f'Getting PubMed IDs (PMID) for articles related to: "{args.search_query}"...')
    search = PubmedSearch(args.email)
    pubmed_ids = search.get_ids(args.search_query, args.start, args.max)
    
    print(f'  {len(pubmed_ids)} PMIDs where found.')
    logging.info(f'{len(pubmed_ids)} PMIDs where found.')

    print('Saving PMIDs to the database...')
    for id in pubmed_ids:
        try:
            database.insert_id(id)
        except sqlite3.IntegrityError:  # will raise on update/retry because of UNIQUE pmid constraint
            continue
    
    pubmed_ids = database.select_new_ids()
    print(f'  {len(pubmed_ids)} PMIDs are marked as NEW in the database.')
    if args.retry:
        failed_ids = database.select_failed_ids()
        pubmed_ids = pubmed_ids + failed_ids
        print(f'  {len(failed_ids)} PMIDs are marked as FAILED in the database.')

    count_success = 0
    count_fail = 0
    count_total = len(pubmed_ids)

    print('Fetching publications...')
    for id in pubmed_ids:
        try:
            publication = search.get_publication(id)
            database.insert_publication(id, publication)
            database.update_fetched_publication(id)
            count_success = count_success + 1
        except (TimeoutError, RuntimeError, TypeError):
            logging.exception(f'ID: {id}')
            database.update_fetched_publication(id, failed=1)
            count_fail = count_fail + 1
        finally:
            print(f'  Successful: {count_success} | Failed: {count_fail} | Remaining: {count_total-count_success-count_fail}', end='\r', flush=True)

    database.connection.close()

    print("\nFinished!")
    logging.info(f'Successful: {count_success} | Failed: {count_fail}')
