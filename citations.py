import re, config, json
from Bio import Entrez
from pprint import pprint

Entrez.email = config.email

# storing sample file as a string
eFetch_sample_JSON = ''
with open('eFetch_sample.json') as data_file:    
    eFetch_sample_JSON = json.load(data_file)

fetch_results = []
webenv = ''
query_key = ''

def search(query):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='20',
                            retmode='xml', 
                            term=query,
                            usehistory='y')
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
	ids = ','.join(id_list)
	handle = Entrez.efetch(db='pubmed',
	                       retmode='xml',
	                       id=ids,
                           webenv=webenv,
                           query_key=query_key)
	results = Entrez.read(handle)
	return results	

# I want to start logging how many papers are found, how many unique addresses,
# how many hits on google places, how many hits on MetaMap
def log():


def first_abstract(results):
    file 

    for r in results:


    print(results[0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0])

# returns a set of strings, each with a different address
# the set of strings with alphanumeric chars only prevents duplication of
# addresses with minor punctuation differences
def unique_addresses(author_list):
    alphanumeric_addresses = set() # addresses with alphanumeric chars only
    result = list()

    for author in author_list:
        for place in author['AffiliationInfo']:
            # splits multiple addresses in a single string into
            # a list of individual addresses
            individual_addresses = place['Affiliation'].split(';')
            for f in individual_addresses:
                formatted_address = _format_address(f)
                # exclude all non-alphanumeric chars, upper case
                alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                if alphanumeric != '' and alphanumeric not in alphanumeric_addresses:
                    result.append(formatted_address)
                    print("New address: " + formatted_address)
                    alphanumeric_addresses.add(alphanumeric)
    print("List of unique addresses: " + str(result))
    return result

# Email addresses are removed to improve success of geocoding using Google Places API
def _remove_email(address):
    # regex for removing nonwhitespace@[alphanum-.]+
    result = re.sub('[\S]+[@][\w.-]+', '', address) 
    reg = re.compile('(email|address|electronic)', re.IGNORECASE)
    result = reg.sub('', result)
    # TODO: make sure comma is also deleted otherwise last thing will be just crap
    return result

# Removes first address lines for addresses over 4 parts long (comma separated)
def _format_address(address):
    without_email = _remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        address_lines = address_lines[-2:] # last 2 elements of the list
    result = (','.join(address_lines))
    return result

###############################
# other ways to structure data!
# institution with all authors (out of scope?)
# paper AND author AND year
# slider for year range then markers appear for author according to search query?
# could be topic (e.g. where/who was writing about DIPG in the last year?)
# could be person (where did papers by Alice Gutteridge originate from in 2014?)
###############################

# separate function that can be called 
def _start_search(query):
    search_results = search(query)
    id_list = search_results['IdList']
    webenv = search_results['WebEnv'] # ID for session
    query_key = search_results['QueryKey'] # ID for query within session

    if not id_list: # if no papers are found
        print("no results!")
        return []
    else:
        papers = fetch_details(id_list)
        # make results global for testing purposes
        # fetch_results = papers
        return papers

# calls BioPython functions esearch and efetch
def get_addresses(query):
    papers = _start_search(query)
    result_list = []

    if not papers: # if no papers are found
        # prompt ask for another query
        print("no results!")
    else:
        print("number of papers: " + str(len(papers)))
        for paper in papers:
            author_list = paper['MedlineCitation']['Article']['AuthorList']
            result_list = result_list + (unique_addresses(author_list))

        first_abstract(result_list)
    return result_list


# write function to make python dict (and therefore JSON) 
# leaner, easier to understand?
# def reorganise(dict):


# testing..
if __name__ == "__main__":
    # _start_search('her2 digital pcr breast cancer')
    first_abstract(eFetch_sample_JSON)
