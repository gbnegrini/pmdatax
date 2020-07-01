import pytest
from pubmed_search import Database, PubmedSearch
from tests.test_pubmed_search import TestPub
import sqlite3
import os

@pytest.fixture
def database_name():
    return 'test'

@pytest.fixture
def database(database_name):
    return Database(database_name)

@pytest.fixture
def pubmed_search():
    return PubmedSearch('gixos30110@lowdh.com') #fake temp mail

def test_db_instance(database_name, database):
    assert os.path.isfile(f'{database_name}.db')

def test_create_connection(database):
    assert type(database._connection) is sqlite3.Connection

def test_create_tables(database):
    cursor = database._connection.cursor()

    cursor.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='pmids'""")
    assert cursor.fetchone()[0] == 1

    cursor.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='publications'""")
    assert cursor.fetchone()[0] == 1

def test_insert_id(database):
    id = 22331878
    database.insert_id(id)
    cursor = database._connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert len(result) == 1
    assert result[0][0] == id
    assert result[0][1] == 0

def test_insert_id_is_unique(database):
    with pytest.raises(sqlite3.IntegrityError):
        id = 22331878
        database.insert_id(id)

def test_insert_publication(database, pubmed_search):
    id = 22331878
    publication = pubmed_search.get_publication(id)
    database.insert_publication(id, publication)

    cursor = database._connection.cursor()
    cursor.execute("""SELECT * FROM publications WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()

    assert len(result) == 1
    assert id == result[0][0]
    assert publication.title == result[0][1]
    assert publication.authors == result[0][2]
    assert publication.journal == result[0][3]
    assert publication.year == str(result[0][4])
    assert publication.month == result[0][5]
    assert publication.day == str(result[0][6])
    assert publication.abstract == result[0][7]

def test_update_fetched_publication(database):
    id = 22331878

    cursor = database._connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert result[0][1] == 0

    database.update_fetched_publication(id)

    cursor = database._connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert result[0][1] == 1
