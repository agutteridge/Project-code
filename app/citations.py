from Bio import Entrez

from app import cache, config

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
    webenv = results['WebEnv'] # ID for session
    query_key = results['QueryKey'] # ID for query within session

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

# Formatting according to # of authors
def get_authors(authors):
    author_list = []
    if len(authors) > 5:
        author_text = authors[0]['LastName'] + ' et al.'
    else:
        author_list = []

        for a in authors:
            # Can be a collective, not a person
            if 'LastName' in a:
                author_list.append(a['LastName'])
            
            if 'CollectiveName' in a:
                author_list.append(a['CollectiveName'])

        author_text = (', ').join(author_list[:-1])
        author_text += ' & ' + author_list[-1]

    return author_text

# Returns either year of publishing or formatted date
def get_date(journal):
    # Dict containing year, month and day 
    if 'PubDate' in journal:
        return journal['PubDate']
    elif 'MedlineDate' in journal:
        return journal['MedlineDate']

# Formatting paper information for frontend as JSON
def format_papers(all_papers):
    results = []

    for p in all_papers:
        author_text = get_authors(p['MedlineCitation']['Article']['AuthorList'])
        date = get_date(p['MedlineCitation']['Article']['Journal']['JournalIssue'])

        result_dict = {'PMID': p['MedlineCitation']['PMID'],
            'title': p['MedlineCitation']['Article']['ArticleTitle'],
            'authors': author_text,
            'date': date,
            'journal': p['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        }

        results.append(result_dict)

    return results

# Returns two lists of strings: the first from the cache with all information,
# and the second from PubMed. 
# Each string is JSON-formatted, and contains data for one paper.
def start_search(query):
    id_list = search(query)['IdList']

    if id_list:
        papers = cache.retrieve(id_list)
        fetch_results = []

        # papers only need to be fetched if not in cache
        if papers['new']:
            fetch_results = fetch_details(papers['new'])

        formatted = format_papers(papers['old'] + fetch_results)

        return {
            'docs': papers['old'],
            'results': fetch_results,
            'formatted': formatted
        }
    else:
        return []
