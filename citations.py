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
    print(author_list[0])

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
                    result.add(formatted_address)
                    alphanumeric_addresses.add(alphanumeric)
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

if __name__ == "__main__":
    print(get_addresses('glioblastoma stem cells'))

