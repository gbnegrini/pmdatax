import argparse
import sqlite3
import os

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

if __name__ == 'main':
    pass