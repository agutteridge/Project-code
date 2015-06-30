import re
from Bio import Entrez
Entrez.email = 'alicegutteridge@gmail.com'

search_term = '"breast cancer"[title] AND her2[title] AND "digital PCR"'

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
            formatted_place = place['Affiliation'].split(';')
            for f in formatted_place:
                # exclude all non-alphanumeric chars
                alphanumeric = re.sub('[\W]', '', f) 
                if alphanumeric not in alphanumeric_addresses:
                    result.add(f)
                    alphanumeric_addresses.add(alphanumeric)
    return result

# some addresses contain emails, these are removed
# before used for Google Places API
def remove_email(string):
    return re.sub('[\s]+[\S]+[@][\w.-]+', '', string)

###############################
# other ways to structure data!
# institution with all authors (out of scope?)
# paper AND author AND year
# slider for year range then markers appear for author according to search query?
# could be topic (e.g. where/who was writing about DIPG in the last year?)
# could be person (where did papers by Alice Gutteridge originate from in 2014?)
###############################

# main method
# calls BioPython functions esearch and efetch
if __name__ == '__main__':
    results = search(search_term)
    id_list = results['IdList']
    webenv = results['WebEnv'] # ID for session
    query_key = results['QueryKey'] # ID for query within session

    if not id_list: # if no results are returned
        print("no results!")
    else:
        papers = fetch_details(id_list)
        # creates tuples of papers with 1-indexed ints 
        for i, paper in enumerate(papers, start=1):
            author_list = paper['MedlineCitation']['Article']['AuthorList']
            print("%d) " % i, end="")
            print(unique_addresses(author_list))