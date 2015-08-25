import re
import json

from Bio import Entrez
import pymongo

from main import db
import config

# Entrez setup
Entrez.email = config.email
webenv = ''
query_key = ''

# Runs a PubMed search on the provided search string
def search(query):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='20',
                            retmode='xml', 
                            term=query,
                            usehistory='y')
    results = Entrez.read(handle)
    return results

# Retrieves citation data for each PubMed ID in a list
def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids,
                           webenv=webenv,
                           query_key=query_key)
    results = Entrez.read(handle)

    return results  

# Formatting paper information for frontend as JSON
def format_papers(all_papers):
    results = []

    for p in all_papers:
        # Formatting according to # of authors
        author_text = ''
        authors = p['MedlineCitation']['Article']['AuthorList']

        if len(authors) > 5:
            author_text = authors[0]['LastName'] + ' et al.'
        else:
            author_list = []
            for a in authors:
                author_list.append(a['LastName'])

            author_text = (', ').join(author_list[:-1])
            author_text += ' & ' + author_list[-1]

        results.append({'PMID': p['MedlineCitation']['PMID'],
            'title': p['MedlineCitation']['Article']['ArticleTitle'],
            'authors': author_text,
            'date': (str(p['MedlineCitation']['Article']['ArticleDate'][0]['Day']) + '/' +
                    str(p['MedlineCitation']['Article']['ArticleDate'][0]['Month']) + '/' +
                    str(p['MedlineCitation']['Article']['ArticleDate'][0]['Year'])),
            'journal': p['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        })

    return results

# Returns two lists of strings: the first from the cache with all information,
# and the second from PubMed. 
# Each string is JSON-formatted, and contains data for one paper.
def start_search(query):
    search_results = search(query)
    id_list = search_results['IdList']
    webenv = search_results['WebEnv'] # ID for session
    query_key = search_results['QueryKey'] # ID for query within session

    if not id_list: # if no papers are found
        print("no results!")
        return ([], [])
    else:
        new_ids = []
        cached_docs = []

        for i in id_list:
            # search for PMID in pubmeddata mongoDB collection
            cursor = db.pubmeddata.find_one( {'MedlineCitation.PMID': i} )

            if not cursor: # if paper not already cached
                new_ids.append(i)
            else:
                cached_docs.append(cursor)

        fetch_results = []

        # papers only need to be fetched if not in cache
        if new_ids:
            fetch_results = fetch_details(new_ids)

        return (cached_docs, fetch_results, format_papers(cached_docs + fetch_results))
