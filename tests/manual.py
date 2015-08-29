import os
import json
import sys
import math

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

def search_100(query):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='100',
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

def getDistanceFromLatLonInKm(lat1,lon1,lat2,lon2):
    r = 6371 # Radius of the earth in km
    dLat = deg2rad(lat2 - lat1)
    dLon = deg2rad(lon2 - lon1)
    a = (math.sin(dLat / 2) * math.sin(dLat / 2) +
        math.cos(deg2rad(lat1)) * math.cos(deg2rad(lat2)) * 
        math.sin(dLon / 2) * math.sin(dLon / 2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c # Distance in km
    return d

def deg2rad(deg):
  return deg * (math.pi / 180)

def format_address_0(address):
    return address

def format_address_1(address):
    return geocode.remove_email(address)

def format_address_2(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    result = (', '.join(without_department))
    return result

def format_address_3(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    result = without_department[0]
    return result

# check > 1
def format_address_4(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    return without_department[1]

# check > 1
def format_address_5(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    return ', '.join(without_department[0:2])

# check > 1
def format_address_6(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    return ', '.join([without_department[0], without_department[-1]])        

# check > 1
def format_address_7(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    return ', '.join([without_department[1], without_department[-1]])

def check(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)

    if len(without_department) > 1:
        True
    else:
        False    

def geocode_run(results):
    with open(os.path.join('./tests/resources', 'format_address_6_results.txt'), 'a') as datafile:
        input_num = 0
        output_num = 0

        for paper in results:
            pmid = str(paper['MedlineCitation']['PMID'])
            datafile.write('\nPMID: ' + pmid + '\n')

            author_list = paper['MedlineCitation']['Article']['AuthorList']

            for author in author_list:
                if 'LastName' in author:
                    datafile.write('\tAuthor: ' + author['LastName'] + '\n')

                elif 'CollectiveName' in author:
                    datafile.write('\tCollective: ' + author['CollectiveName'] + '\n')

                for place in author['AffiliationInfo']:
                    individual_addresses = place['Affiliation'].split(';')

                    for f in individual_addresses:
                        if (check(f)):
                            datafile.write('\t\tInput address: ' + f + '\n')
                            input_num += 1
                            formatted_address = format_address_6(f) # change format_address
                            place = geocode.get_location(formatted_address)
                            datafile.write('\t\tFormatted address: ' + formatted_address + '\n')                            
        
                            if place:
                                datafile.write('\t\tOutput address: ' + place['name'] + '\n')
                                output_num += 1
                            else:
                                datafile.write('\t\tOutput address: NONE' + '\n')

        datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%\n')
    datafile.close()

# Writes first 20 papers of chosen universities
def init_spec(unis):
    with open(os.path.join('./tests/resources', 'eFetch_unis.json'), 'a') as datafile:
        print('opened')
        json_list = []

        for u in unis:
            id_list = citations._search(u)['IdList']
            fetch_results = citations.fetch_details(id_list)

            for fr in fetch_results:
                json_list.append(fr)

        datafile.write(json.dumps(json_list))
    
    datafile.close()

# Writes first 100 papers from PubMed retrieved with query 'biology' to JSON file 
def init_sens():
    with open(os.path.join('./tests/resources', 'eFetch_biology.json'), 'a') as datafile:
        print('opened')

        id_list = search_100('biology')['IdList']
        fetch_results = citations.fetch_details(id_list)
        json_list = []

        for fr in fetch_results:
            json_list.append(fr)

        datafile.write(json.dumps(json_list))
    datafile.close()

if __name__ == "__main__":
    init_spec([
        'Massachusetts Institute of Technology[Affiliation]',
        'University College London[Affiliation]',
        'University of Cape Town[Affiliation]',
        'École Normale Supérieure[Affiliation]',
        'University of Tokyo[Affiliation]'
    ])
    # geocode_run(fake_json('eFetch_biology.json'))
    # print('no method chosen')
