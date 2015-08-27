import os
import json
import sys

from Bio import Entrez

# configuring relative imports
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from app import citations, geocode, cache, config

# Entrez setup
Entrez.email = config.email
webenv = ''
query_key = ''

def citations_search(query):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='101',
                            retmode='xml', 
                            term=query,
                            usehistory='y')
    results = Entrez.read(handle)
    webenv = results['WebEnv'] # ID for session
    query_key = results['QueryKey'] # ID for query within session

    return results

def fake_json(filename):
    with open(os.path.join('./tests/resources', filename), 'r') as datafile:
        obj = json.load(datafile)
        datafile.close()
        return obj

def format_address_0(address):
    return address

def format_address_1(address):
    return geocode.remove_email(address)

def format_address_2(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    result = (', '.join(without_department))
    return result

def geocode_run(results):
    with open(os.path.join('./tests/resources', 'format_address_2_results.txt'), 'a') as datafile:
        input_num = 0
        output_num = 0

        for paper in results:
            pmid = str(paper['MedlineCitation']['PMID'])
            datafile.write('\nPMID: ' + pmid + '\n')

            if 'AuthorList' in paper['MedlineCitation']['Article']: 
                author_list = paper['MedlineCitation']['Article']['AuthorList']

                for author in author_list:
                    if 'LastName' in author:
                        datafile.write('\tAuthor: ' + author['LastName'] + '\n')

                    if 'CollectiveName' in author:
                        datafile.write('\tCollective: ' + author['CollectiveName'] + '\n')

                    for place in author['AffiliationInfo']:
                        individual_addresses = place['Affiliation'].split(';')

                        for f in individual_addresses:
                            datafile.write('\t\tInput address: ' + f + '\n')
                            input_num += 1
                            formatted_address = format_address_2(f) # change format_address
                            place = geocode.get_location(formatted_address)
                            datafile.write('\t\tFormatted address: ' + formatted_address + '\n')                            
        
                            if place:
                                datafile.write('\t\tOutput address: ' + place['name'] + '\n')
                                output_num += 1
                            else:
                                datafile.write('\t\tOutput address: NONE' + '\n')
            else:
                datafile.write('\tNO AUTHORS LISTED.\n')

        datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%')
    datafile.close()


# Writes first 100 papers from PubMed retrieved with query 'biology' to JSON file 
def init():
    with open(os.path.join('./tests/resources', 'eFetch_biology.json'), 'a') as datafile:
        print('opened')

        id_list = citations_search('biology')['IdList']
        fetch_results = citations.fetch_details(id_list)
        json_list = []

        for fr in fetch_results:
            if 'MedlineCitation' in fr:
                json_list.append(fr)

        datafile.write(json.dumps(json_list))
    datafile.close()

if __name__ == "__main__":
    # count(fake_json('eFetch_biology.json'))
    init()
    # geocode_run(fake_json('eFetch_biology.json'))
    # print('no method chosen')
