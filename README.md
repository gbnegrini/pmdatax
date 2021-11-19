
<img src="https://i.ibb.co/df5ynrX/undraw-text-files-au1q.png" alt="undraw-text-files-au1q" border="0">

# PMdataX - PubMed data eXtraction
> A command-line tool to fetch data from articles available on PubMed

[![Build Status](https://travis-ci.com/gbnegrini/pmdatax.svg?token=QivxGnspzs4tLzQsGSae&branch=master)](https://travis-ci.com/gbnegrini/pmdatax)
[![Generic badge](https://img.shields.io/badge/python-3.6|3.7|3.8-blue.svg)](https://shields.io/)
![LICENSE](https://img.shields.io/github/license/gbnegrini/pmdatax)

For a given search query, PMdataX uses [Bio.Entrez](https://biopython.org/docs/1.74/api/Bio.Entrez.html) to get the PubMed IDs from related articles. Each retrieved PubMed ID is then used to fetch its publication data with [pubmed-lookup](https://github.com/mfcovington/pubmed-lookup). All data is saved to a SQLite database and can be easily exported to a csv file. Currently, these data include: <b>PMID, title, authors, journal, year, month, day, abstract, publication type, DOI and URL</b>.

## Requirements
- python>=3.7
- biopython==1.76
- pubmed-lookup==0.2.3
- pandas==1.3.2
- pytest==6.2.5 (optional, required only for code testing)

## Installation
Clone this repository, make sure you have installed the requirements mentioned above and you are good to go.
This can be easily achieved by running: 
```
git clone https://github.com/gbnegrini/pmdatax.git
cd pmdatax/
pip install -r requirements.txt
```

Alternatively, a **[Docker image](Dockerfile)** is available.

## Quick start
Let's fetch some data about COVID-19 publications. Note we're using the `--max` flag in this example so it doesn't take long to run.

```bash
python pmdatax.py fetch "covid-19 OR sars-cov2" my_covid_database your_email@somethingemail.com --max 10
```

If there are any failed publications, try again with: 

```bash
python pmdatax.py retry my_covid_database your_email@somethingemail.com
```

Once everything is finished, export to csv if you don't want to work with SQL:

```bash
python pmdatax.py export my_covid_database
```

### Using Docker

You can use the image to run the application inside a container.

```bash
docker build -t pmdatax .

docker run --rm -it -v $(pwd):/app pmdatax python pmdatax.py fetch "covid-19 OR sars-cov2" my_covid_database your_email@somethingemail.com --max 10
```

## Detailed usage
PMdataX supports three different commands:
- fetch
- retry
- export
```
python pmdatax.py -h

usage: pmdatax.py [-h] {fetch,retry,export} ...

PMdataX - Fetch data from articles available on PubMed!

positional arguments:
  {fetch,retry,export}
    fetch               Fetch PubMed data
    retry               Try to fetch again any PMIDs marked as failed in the database
    export              Export publication data from the database to a <database_name>.csv file

optional arguments:
  -h, --help            show this help message and exit
````
### fetch
`fetch` is the primary command and it is used to fetch data from articles related to a given search query.

This operation can take several hours depending on how many articles are available. You can kill the process any time you want and run the same command to fetch the remaining articles data. Running the same command, <i>i.e.</i> same search query and database name, is also a way to keep your records up to date.

```
python pmdatax.py fetch -h

usage: pmdatax.py fetch [-h] [-s START] [-m MAX] search_query database_name email

positional arguments:
  search_query          Search query for PubMed articles. Use quotation marks to write more than a single word
  database_name         Create or connect to an existing <database_name>.db file
  email                 Identify yourself so NCBI can contact you in case of excessive requests

optional arguments:
  -h, --help            show this help message and exit
  -s START, --start START
                        Sequential index of the first PMID from retrieved set (default: 0)
  -m MAX, --max MAX     Total number of PMIDs from the retrieved set. Max: 100,000 records (default: 100000)
```
### retry
`retry` is used to try to fetch again any data that failed to be retrieved.

Sometimes an error occurs when fetching data. Errors are mainly caused by Internet connection instabillity or missing fields from source. When this happens, PMdataX marks the PMID as failed in the database and logs the error in a <database_name>_error.log file. `retry` will try to fetch only these articles, but keep in mind that if the error was due to missing or incorrect data from source then it will fail again.

```
python pmdatax.py retry -h

usage: pmdatax.py retry [-h] database_name email

positional arguments:
  database_name  Connect to an existing <database_name>.db file
  email          Identify yourself so NCBI can contact you in case of excessive requests

optional arguments:
  -h, --help     show this help message and exit
```

### export
`export` command will easily export the data to a csv file.

```
python pmdatax.py export -h

usage: pmdatax.py export [-h] database_name

positional arguments:
  database_name  Connect to an existing <database_name>.db file

optional arguments:
  -h, --help     show this help message and exit
```

## Output data
PMdataX will save the extracted data to a SQLite database.
The database can be queried using standard SQL sintax. For more information, please consult the [SQLite page](https://www.sqlite.org/lang.html). Moreover, you can export the publications data to a csv file.

### SQLite
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
    - type TEXT
    - doi TEXT
    - url TEXT

### CSV
csv file will look like this:

| pmid | title | authors | journal | year | month | day | abstract | type | doi | url |
| --------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- | --------- |
| 22331878 | Arabidopsis synchronizes jasmonate-mediated defense with insect circadian behavior. | "Goodspeed D, Chehab EW, Min-Venditti A, Braam J, Covington MF" | Proc Natl Acad Sci U S Arnal | 2012 | 3 | 20 | "Diverse life forms have evolved internal clocks[...]" | Journal Article | 10.1073/pnas.1116368109 | https://www.pnas.org/content/109/12/4674 |