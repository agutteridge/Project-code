import re
from Bio import Entrez
Entrez.email = 'alicegutteridge@gmail.com'

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
	                       id=ids
                           #webenv and query cause NameError?
                           )
	results = Entrez.read(handle)
	return results	

# returns a set of strings, each with a different address
# the set of strings with alphanumeric chars only prevents duplication of
# addresses with minor punctuation differences
def unique_addresses(author_list):
    alphanumeric_addresses = set() # addresses with alphanumeric chars only
    result = set()

    for author in author_list:
        for place in author['AffiliationInfo']:
            # splits multiple addresses in a single string into
            # a list of individual addresses
            individual_addresses = place['Affiliation'].split(';')
            for f in individual_addresses:
                # exclude all non-alphanumeric chars
                alphanumeric = re.sub('[\W]', '', f) 
                if alphanumeric not in alphanumeric_addresses:
                    formatted_address = _format_address(f)
                    result.add(formatted_address)
                    alphanumeric_addresses.add(alphanumeric)
    return result

# some addresses contain emails, these are removed
# before used for Google Places API
def _remove_email(address):
    # regex for removing nonwhitespace@[alphanum-.]+
    result = re.sub('[\S]+[@][\w.-]+', '', address)
    word_list =  ['email', 'address', 'electronic']
    insensitive_hippo = re.compile(re.escape('hippo'), re.IGNORECASE)
    for w in word_list:
        result = result.replace(w, '')
    # make sure comma is also deleted otherwise last thing will be just crap
    return result

# removes first address lines for addresses over 4 parts long (comma separated)
def _format_address(address):
    address_lines = address.split(',')
    if len(address_lines) > 4:
        address_lines = address_lines[-4:] # last 4 elements of the list
    result = _remove_email(','.join(address_lines))
    return result

###############################
# other ways to structure data!
# institution with all authors (out of scope?)
# paper AND author AND year
# slider for year range then markers appear for author according to search query?
# could be topic (e.g. where/who was writing about DIPG in the last year?)
# could be person (where did papers by Alice Gutteridge originate from in 2014?)
###############################

# calls BioPython functions esearch and efetch
def get_addresses(query):
    results = search(query)
    id_list = results['IdList']
    webenv = results['WebEnv'] # ID for session
    query_key = results['QueryKey'] # ID for query within session
    result_set = set()

    if not id_list: # if no results are returned
        print("no results!")
    else:
        papers = fetch_details(id_list)
        # creates tuples of papers with 1-indexed ints
        for i, paper in enumerate(papers, start=1):
            author_list = paper['MedlineCitation']['Article']['AuthorList']
            result_set = result_set.union(unique_addresses(author_list))

    return result_set        