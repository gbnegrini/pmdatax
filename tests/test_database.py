import pytest
from pubmed_search import Database
import sqlite3
import os

@pytest.fixture
def database_name():
    return 'test'

@pytest.fixture
def example_database(database_name):
    return Database(database_name)

def test_db_instance(database_name, example_database):
    assert os.path.isfile(f'{database_name}.db')

def test_create_connection(example_database):
    assert type(example_database._connection) is sqlite3.Connection

def test_create_tables(example_database):
    cursor = example_database._connection.cursor()

    cursor.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='pmids'""")
    assert cursor.fetchone() == 1

    cursor.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='publications'""")
    assert cursor.fetchone() == 1