import os
import json
import sys
import math
import re

from Bio import Entrez

# configuring relative imports
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from app import citations, geocode, cache, config
import graphs

# Entrez setup
Entrez.email = config.email
webenv = ''
query_key = ''

def search_num(query, num):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax=num,
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

# check > 0
def format_address_0(address):
    return address

# check > 0
def format_address_1(address):
    return geocode.remove_email(address)

# check > 0
def format_address_2(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    lines_list.append(len(without_department))
    result = (', '.join(without_department))
    return result

# check > 0
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
    address_lines = address.split(',')
    if len(address_lines) > 0:
        return True
    else:
        return False    

# def geocode_specificity(results):
#     coords = fake_json('spec_results.json')
#     coordinate_list = []

#         for paper in results:
#             pmid = str(paper['MedlineCitation']['PMID'])
#             datafile.write('\nPMID: ' + pmid + '\n')

#             author_list = paper['MedlineCitation']['Article']['AuthorList']
#             alphanumeric_addresses = set() # do not analyse duplicate addresses per paper

#             for author in author_list:
#                 if 'LastName' in author:
#                     datafile.write('\tAuthor: ' + author['LastName'] + '\n')

#                 elif 'CollectiveName' in author:
#                     datafile.write('\tCollective: ' + author['CollectiveName'] + '\n')

#                 for place in author['AffiliationInfo']:
#                     individual_addresses = place['Affiliation'].split(';')

#                     coords_number = 0
#                     for f in individual_addresses:
#                         formatted_address = format_address_1(f) # change format_address
#                         alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                        
#                         if alphanumeric != '' and alphanumeric not in alphanumeric_addresses and check(f):
#                             place = geocode.get_location(formatted_address)
#                             datafile.write('\t\tFormatted address: ' + formatted_address + '\n')                            
        
#                             if place:
#                                 datafile.write('\t\tOutput address: ' + place['name'] + '\n')
#                                 this_coord = place['geometry']['location']
#                                 coordinate_list.append()
#                                 distance = getDistanceFromLatLonInKm(this_coord['lat'],
#                                     this_coord['long'],
#                                     coords[coords_number]['lat'],
#                                     coords[coords_number]['long'])
#                                 datafile.write('\t\tDistance from target: ' + str(distance) + ' km\n')
#                                 output_num += 1
#                             else:
#                                 coordinate_list.append({'coordinates': False})
#                                 datafile.write('\t\tOutput address: NONE' + '\n')

#                         coords_number += 1 # move to next coord in list regardless

#         if input_num > 0:
#             datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%\n')
#     datafile.close()

#     with open(os.path.join('./tests/resources', 'format_address_01_coords.json'), 'a') as datafile:
#         datafile.write(json.dumps(coordinate_list))
#     datafile.close()

lines_list = []
def geocode_sensitivity(results):
    with open(os.path.join('./tests/resources', '02_sensitivity.txt'), 'a') as datafile:
        input_num = 0
        output_num = 0

        for paper in results:
            pmid = str(paper['MedlineCitation']['PMID'])
            datafile.write('\nPMID: ' + pmid + '\n')

            author_list = paper['MedlineCitation']['Article']['AuthorList']
            alphanumeric_addresses = set() # do not analyse duplicate addresses per paper

            for author in author_list:
                for place in author['AffiliationInfo']:
                    individual_addresses = place['Affiliation'].split(';')

                    for f in individual_addresses:
                        formatted_address = format_address_2(f) # change format_address
                        alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                        
                        if alphanumeric != '' and alphanumeric not in alphanumeric_addresses and check(f):
                            datafile.write('\tInput address: ' + f + '\n')
                            input_num += 1
                            place = geocode.get_location(formatted_address)
                            datafile.write('\tFormatted address: ' + formatted_address + '\n')                            
        
                            if place:
                                datafile.write('\tOutput address: ' + place['name'] + '\n')
                                output_num += 1
                            else:
                                datafile.write('\tOutput address: NONE' + '\n')

        if input_num > 0:
            datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%\n' + str(countries_dict))
    datafile.close()
    graphs.plot(lines_list)

def init_specificity():
    id_list = (search_num('USA[Affiliation]', 1)['IdList'] + 
        search_num('China[Affiliation]', 1)['IdList'] + 
        search_num('UK[Affiliation]', 1)['IdList'] + 
        search_num('Japan[Affiliation]', 1)['IdList'] +
        search_num('Germany[Affiliation]', 1)['IdList'] + 
        search_num('Italy[Affiliation]', 1)['IdList'] + 
        search_num('Canada[Affiliation]', 1)['IdList'] + 
        search_num('France[Affiliation]', 1)['IdList'] + 
        search_num('Australia[Affiliation]', 1)['IdList'] + 
        search_num('India[Affiliation]', 1)['IdList'])

    unique_ids = []
    for i in id_list:
        if i not in unique_ids:
            unique_ids.append(i)

    fetch_results = citations.fetch_details(unique_ids)
    json_list = []

    for fr in fetch_results:
        json_list.append(fr)

    with open(os.path.join('./tests/resources', 'eFetch_specificity.json'), 'w') as datafile:
        print('opened, writing...')
        datafile.write(json.dumps(json_list))
    datafile.close()

# Retrieves unique papers according to affiliate countries,
# according to % share of citations
def init_sensitivity():
    id_list = (search_num('USA[Affiliation]', 14)['IdList'] + 
        search_num('China[Affiliation]', 4)['IdList'] + 
        search_num('UK[Affiliation]', 4)['IdList'] + 
        search_num('Japan[Affiliation]', 3)['IdList'] +
        search_num('Germany[Affiliation]', 3)['IdList'] + 
        search_num('Italy[Affiliation]', 2)['IdList'] + 
        search_num('Canada[Affiliation]', 2)['IdList'] + 
        search_num('France[Affiliation]', 2)['IdList'] + 
        search_num('Australia[Affiliation]', 1)['IdList'] + 
        search_num('India[Affiliation]', 1)['IdList'])

    unique_ids = []
    for i in id_list:
        if i not in unique_ids:
            unique_ids.append(i)

    fetch_results = citations.fetch_details(unique_ids)

    with open(os.path.join('./tests/resources', 'eFetch_sensitivity.json'), 'w') as datafile:
        print('opened, writing...')
        datafile.write(json.dumps(fetch_results))
    datafile.close()

if __name__ == "__main__":
    # init_sensitivity()
    # init_specificity()
    geocode_sensitivity(fake_json('eFetch_sensitivity.json'))
    # print('no method chosen')
