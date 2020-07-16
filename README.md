
<img src="https://i.ibb.co/df5ynrX/undraw-text-files-au1q.png" alt="undraw-text-files-au1q" border="0">

# PubMed data extraction
> A command-line script to fetch data from PubMed articles

[![Build Status](https://travis-ci.com/gbnegrini/pubmed-text-mining.svg?token=QivxGnspzs4tLzQsGSae&branch=master)](https://travis-ci.com/gbnegrini/pubmed-text-mining)
[![Generic badge](https://img.shields.io/badge/python-3.6|3.7|3.8-blue.svg)](https://shields.io/)
![LICENSE](https://img.shields.io/github/license/gbnegrini/pubmed-text-mining)

For a given search query, the script uses [Bio.Entrez package](https://biopython.org/docs/1.74/api/Bio.Entrez.html) to get the PubMed IDs from related articles. Each retrieved PubMed ID is then used to fetch its publication data with [pubmed-lookup package](https://github.com/mfcovington/pubmed-lookup). 

## Requirements
- python>=3.6
- biopython==1.76
- pubmed-lookup==0.2.3
- pytest==5.4.3 (optional, required only for code testing)

## Usage
```
python pubmed_extraction.py <email> <search_query> [-h] [-s START] [-m MAX] [-o OUTPUT] [-r]

Get publication data from PubMed.

positional arguments:
  email                 Identify yourself so NCBI can contact you in case of excessive requests
  search_query          Search term for Pubmed database

optional arguments:
  -h, --help            show this help message and exit
  -s START, --start START
                        Sequential index of the first PMID from retrieved set. (default: 0)
  -m MAX, --max MAX     Total number of PMIDs from the retrieved set. Max: 100,000 records. (default: 100000)
  -o OUTPUT, --output OUTPUT
                        Specify an existing database. (default: create or connect to a database using the search query as a file name).    
  -r, --retry           Try to fetch again any PMIDs marked as failed in the database.

````

## Output data
The script will save the extracted data to a SQLite database.
The database can be queried using standard SQL sintax. For more information, please consult the [SQLite page](https://www.sqlite.org/lang.html).

There are two tables organized as follows:
- pmids

    - pmid INTEGER NOT NULL PRIMARY KEY
    - new BOOLEAN NOT NULL CHECK (new IN (0,1))
    - failed BOOLEAN NOT NULL CHECK (failed IN (0,1))

- publications
    - pmid INTEGER NOT NULL
    - title TEXT
    - authors TEXT
    - journal TEXT
    - year INTEGER
    - month INTEGER
    - day INTEGER
    - abstract TEXT
