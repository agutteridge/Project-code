import os
import json
import sys
import math
import re
import statistics

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
    d = r * c # Distance in km
    return d

def deg2rad(deg):
  return deg * (math.pi / 180)

def format_address_00(address):
    return address

def format_address_01(address):
    return geocode.remove_email(address)

def format_address_02(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    result = (', '.join(without_department))
    return result

def format_address_03(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 0:
        return without_department[0]
    else:
        return ''

def format_address_04(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 1:
        return without_department[1]
    else:
        return ''

def format_address_05(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 1:
        return ', '.join(without_department[0:2])
    else:
        return ''

def format_address_06(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 1:
        return ', '.join([without_department[0], without_department[-1]])        
    else:
        return ''

def format_address_07(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 2:
        return ', '.join([without_department[1], without_department[-1]])
    else:
        return ''

def format_address_08(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 2:
        return without_department[2]
    else:
        return ''

def format_address_09(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 3:
        return ', '.join([without_department[2], without_department[-1]])
    else:
        return ''

def format_address_10(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 1:
        return ', '.join(without_department[-2:])
    else:
        return ''

def format_address_11(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 2:
        return ', '.join(without_department[-3:])
    else:
        return ''

def format_address_12(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 2:
        return ', '.join(without_department[:3])
    else:
        return ''

def format_address_13(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 3:
        return ', '.join(without_department[-4:])
    else:
        return ''

def format_address_14(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 3:
        return ', '.join(without_department[:4])
    else:
        return ''

def check(address):
    if len(address) > 1:
        return True
    else:
        return False    

def geocode_specificity(results):
    with open(os.path.join('./tests/resources', '14_specificity.txt'), 'a') as datafile:
        input_num = 0
        output_num = 0
        distance_list = []

        for paper in results:
            pmid = str(paper['MedlineCitation']['PMID'])
            datafile.write('\nPMID: ' + pmid + '\n')

            author_list = paper['MedlineCitation']['Article']['AuthorList']

            for author in author_list:
                for a in author['AffiliationInfo']:
                    
                    if 'lat' in a: # coordinates have been found manually
                        formatted_address = format_address_14(a['Affiliation']) # change format_address
                        if formatted_address:
                            datafile.write('\tInput address: ' + a['Affiliation'] + '\n')
                            place = geocode.get_location(formatted_address)
                            datafile.write('\t\tFormatted address: ' + formatted_address + '\n')                            
        
                            if place:
                                input_num += 1
                                datafile.write('\tOutput address: ' + place['name'] + '\n')
                                datafile.write('\tResult coordinates: ' + str(place['geometry']['location']) + '\n')
                                distance = getDistanceFromLatLonInKm(
                                    place['geometry']['location']['lat'],
                                    place['geometry']['location']['lng'],
                                    a['lat'],
                                    a['lng'])
                                distance_list.append(distance)

                                if distance < 5:
                                    output_num += 1
                                
                                datafile.write('\t\tDistance from target: ' + str(distance) + ' km\n')
                            else:
                                datafile.write('\t\tOutput address: NONE' + '\n')

        datafile.write(str(distance_list) + '\n\n')
        datafile.write('\nMEAN DISTANCE: ' + str(statistics.mean(distance_list)) + ' km\n')
        datafile.write('\n1 STANDARD DEVIATION: ' + str(statistics.stdev(distance_list)) + ' km\n')

        if input_num > 0:
            datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%\n')    
    datafile.close()

def geocode_sensitivity(results):
    with open(os.path.join('./tests/resources', '14_sensitivity.txt'), 'a') as datafile:
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
                        formatted_address = format_address_14(f) # change format_address
                        alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                        
                        if alphanumeric != '' and alphanumeric not in alphanumeric_addresses and formatted_address:
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
            datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%\n')
    datafile.close()

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
    # geocode_sensitivity(fake_json('eFetch_sensitivity.json'))
    # geocode_specificity(fake_json('eFetch_specificity.json'))
    print('no method chosen')
