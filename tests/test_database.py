import pytest
from pmdatax import Database, PubmedSearch
import sqlite3
import os

@pytest.fixture
def database_name():
    return 'test'

@pytest.fixture
def empty_database(database_name):
    db = Database(f'{database_name}.db')
    yield db
    db.connection.close()
    os.remove(f'{database_name}.db')

@pytest.fixture
def pubmed_search():
    return PubmedSearch('gixos30110@lowdh.com') #fake temp mail

@pytest.fixture
def dummy_database(database_name, pubmed_search):
    db = Database(f'{database_name}.db')
    db.insert_id(22331878)
    db.insert_publication(22331878, pubmed_search.get_publication(22331878))
    db.insert_id(22331879)
    db.insert_publication(22331879, pubmed_search.get_publication(22331879))
    yield db
    db.connection.close()
    os.remove(f'{database_name}.db')

def test_connection(empty_database):
    assert type(empty_database.connection) is sqlite3.Connection

def test_create_tables(empty_database):
    cursor = empty_database.connection.cursor()

    cursor.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='pmids'""")
    assert cursor.fetchone()[0] == 1

    cursor.execute("""SELECT count(name) FROM sqlite_master WHERE type='table' AND name='publications'""")
    assert cursor.fetchone()[0] == 1

def test_insert_id(empty_database):
    id = 22331878
    empty_database.insert_id(id)
    cursor = empty_database.connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert len(result) == 1
    assert result[0][0] == id
    assert result[0][1] == 1
    assert result[0][2] == 0

def test_insert_id_is_unique(dummy_database):
    with pytest.raises(sqlite3.IntegrityError):
        dummy_database.insert_id(22331878)

def test_insert_publication(dummy_database, pubmed_search):
    id = 22331880
    dummy_database.insert_id(id)
    publication = pubmed_search.get_publication(id)
    dummy_database.insert_publication(id, publication)

    cursor = dummy_database.connection.cursor()
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

def test_insert_publication_without_foreign_key(empty_database, pubmed_search):
    id = 22331878
    publication = pubmed_search.get_publication(id)
    with pytest.raises(sqlite3.IntegrityError):
        empty_database.insert_publication(id, publication)

def test_update_fetched_publication(dummy_database):
    id = 22331878
    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert len(result) == 1
    assert result[0][0] == id
    assert result[0][1] == 1
    assert result[0][2] == 0

    dummy_database.update_fetched_publication(id)

    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert len(result) == 1
    assert result[0][0] == id
    assert result[0][1] == 0
    assert result[0][2] == 0

def test_update_fetched_publication_failed(dummy_database):
    id = 22331879
    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert len(result) == 1
    assert result[0][0] == id
    assert result[0][1] == 1
    assert result[0][2] == 0

    dummy_database.update_fetched_publication(id, failed=1)

    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT * FROM pmids WHERE pmid = ?""", [id])
    result = cursor.fetchall()
    cursor.close()
    assert len(result) == 1
    assert result[0][0] == id
    assert result[0][1] == 0
    assert result[0][2] == 1

def test_select_new_ids(dummy_database):
    new_pmids = dummy_database.select_new_ids()
    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT pmid FROM pmids WHERE new = 1""")
    result = cursor.fetchall()
    cursor.close()

    assert len(result) == 2
    assert len(new_pmids) == 2
    assert new_pmids[0] == result[0][0]
    assert new_pmids[1] == result[1][0]

    dummy_database.update_fetched_publication(22331879, new=0)

    new_pmids = dummy_database.select_new_ids()
    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT pmid FROM pmids WHERE new = 1""")
    result = cursor.fetchall()
    cursor.close()

    assert len(result) == 1
    assert len(new_pmids) == 1
    assert new_pmids[0] == result[0][0]

def test_select_failed_ids(dummy_database):
    dummy_database.update_fetched_publication(22331879, new=0, failed=1)

    failed_pmids = dummy_database.select_failed_ids()
    cursor = dummy_database.connection.cursor()
    cursor.execute("""SELECT pmid FROM pmids WHERE failed = 1""")
    result = cursor.fetchall()
    cursor.close()

    assert len(result) == 1
    assert len(failed_pmids) == 1
    assert failed_pmids[0] == result[0][0]

def test_export_to_csv(dummy_database, database_name):
    dummy_database.export_to_csv(f'{database_name}.csv')
    assert os.path.exists(f'{database_name}.csv') == True
    #os.remove(f'{database_name}.csv')