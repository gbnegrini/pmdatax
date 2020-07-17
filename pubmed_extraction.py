#!/usr/bin/env python
import argparse
import sqlite3
import os
from Bio import Entrez
from pubmed_lookup import PubMedLookup, Publication
import logging
import sys
import pandas as pd

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

    def export_to_csv(self, output_name: str):
        table = pd.read_sql_query("SELECT * FROM publications", self.connection)
        table.to_csv(output_name, index=None, header=True)

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

def _parse_args(args: list) -> argparse.ArgumentParser.parse_args:
    parser = argparse.ArgumentParser(description='Get publication data from PubMed articles! PMID, title, authors, journal, year, month, day and abstract.')
    subparsers = parser.add_subparsers(dest="command")

    fetch = subparsers.add_parser('fetch', help='Fetch PubMed data')
    fetch.add_argument('search_query', type=str,help='Search query for PubMed articles. Use quotation marks to write more than a single word')
    fetch.add_argument('database_name', type=str, help='Create or connect to an existing <database_name>.db file')
    fetch.add_argument('email', type=str, help='Identify yourself so NCBI can contact you in case of excessive requests')
    fetch.add_argument('-s', '--start', type=int, help='Sequential index of the first PMID from retrieved set (default: 0)', default=0)
    fetch.add_argument('-m', '--max', type=int, help='Total number of PMIDs from the retrieved set. Max: 100,000 records (default: 100000)', default=100000)
    
    retry = subparsers.add_parser('retry', help='Try to fetch again any PMIDs marked as failed in the database')
    retry.add_argument('database_name', type=str,help='Connect to an existing <database_name>.db file')
    retry.add_argument('email', type=str, help='Identify yourself so NCBI can contact you in case of excessive requests')
    
    export = subparsers.add_parser('export', help='Export publication data from the database to a <database_name>.csv file')
    export.add_argument('database_name', type=str,help='Connect to an existing <database_name>.db file')
    
    return parser.parse_args()

if __name__ == '__main__':

    args = _parse_args(sys.argv[1:])
    logging.basicConfig(filename=f'{args.database_name}_error.log',format='%(asctime)s - %(message)s', level=logging.INFO)

    if args.command == 'fetch' or args.command == 'retry':
        
        database = Database(f'{args.database_name}.db')
        search = PubmedSearch(args.email)

        if args.command == 'fetch':

            print(f'Getting PubMed IDs (PMID) for articles related to: "{args.search_query}"...')
            pubmed_ids = search.get_ids(args.search_query, args.start, args.max)
            print(f'  {len(pubmed_ids)} PMIDs where found.')

            print('Saving PMIDs to the database...')
            for id in pubmed_ids:
                try:
                    database.insert_id(id)
                except sqlite3.IntegrityError:  # will raise on update because of UNIQUE pmid constraint
                    continue
            
            pubmed_ids = database.select_new_ids()
            print(f'  {len(pubmed_ids)} PMIDs are marked as NEW in the database.')

        elif args.command == 'retry':
            pubmed_ids = database.select_failed_ids()
            print(f'  {len(pubmed_ids)} PMIDs are marked as FAILED in the database.')

        print('Fetching publications...')
        count_success = 0
        count_fail = 0
        count_total = len(pubmed_ids)
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
            
    elif args.command == 'export':
        database = Database(f'{args.database_name}.db')
        database.export_to_csv(f'{args.database_name}.csv')

    database.connection.close()
    
    print("\nFinished!")
    