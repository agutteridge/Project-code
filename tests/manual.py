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

from app import citations, geocode, cache, config, metamap

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

# No formatting
def format_address_00(address):
    return address

# Email address removed
def format_address_01(address):
    return geocode.remove_email(address)

# Department removed
def format_address_02(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    result = (', '.join(without_department))
    return result

# First line removed
def format_address_03(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 1:
        return ', '.join(address_lines[1:len(address_lines)])
    else:
        return ''

# All but first line, department removed
def format_address_04(address):
    without_email = geocode.remove_email(address)
    without_department = geocode.remove_dept(without_email)
    if len(without_department) > 1:
        return ', '.join(without_department[1:len(without_department)])
    else:
        return ''

# Last 2 lines
def format_address_05(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 1:
        return ', '.join(address_lines[-2:])
    else:
        return ''

# Second and third lines only
def format_address_06(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        return ', '.join(address_lines[1:3])
    else:
        return ''

# Second and last lines
def format_address_07(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        return ', '.join([address_lines[1], address_lines[-1]])
    else:
        return ''

# Last 3 lines
def format_address_08(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        return ', '.join(address_lines[-3:])
    else:
        return ''

# All but first 2 lines
def format_address_09(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        return ', '.join(address_lines[2:len(address_lines)])
    else:
        return ''

def format_address_10(address):
    without_email = geocode.remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        return ', '.join(address_lines[1:len(address_lines)-1])
    else:
        return ''    

def geocode_acc(results):
    with open(os.path.join('./tests/resources', '10_acc.txt'), 'a') as datafile:
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
                        formatted_address = format_address_10(a['Affiliation']) # change format_address
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

def geocode_hit(results):
    with open(os.path.join('./tests/resources', '10_hit.txt'), 'a') as datafile:
        input_num = 0
        output_num = 0
        alphanumeric_addresses = set() # do not analyse duplicate addresses

        for paper in results:
            pmid = str(paper['MedlineCitation']['PMID'])
            datafile.write('\nPMID: ' + pmid + '\n')

            author_list = paper['MedlineCitation']['Article']['AuthorList']

            for author in author_list:
                for place in author['AffiliationInfo']:
                    individual_addresses = place['Affiliation'].split(';')

                    for f in individual_addresses:
                        formatted_address = format_address_10(f) # change format_address
                        alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                        
                        if alphanumeric != '' and alphanumeric not in alphanumeric_addresses and formatted_address:
                            alphanumeric_addresses.add(alphanumeric)
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
def init_hit():
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

    with open(os.path.join('./tests/resources', 'eFetch_hit.json'), 'w') as datafile:
        print('opened, writing...')
        datafile.write(json.dumps(fetch_results))
    datafile.close()


def doublecheck():
    place = """TN, USA"""
    print(geocode.get_location(place))

def acc_tiered():
    distance_list = []

    distance_list.append(getDistanceFromLatLonInKm(53.5890864, 9.9759963,
                                                   53.596819, 10.153465))
    distance_list.append(getDistanceFromLatLonInKm(35.1555334, -90.0455115,
                                                   35.155527, -90.045487))
    distance_list.append(getDistanceFromLatLonInKm(59.3481484, 18.0236578,
                                                   59.346709, 18.039748))
    distance_list.append(getDistanceFromLatLonInKm(40.022621, 116.2860293,
                                                   40.022621, 116.286029))
    distance_list.append(getDistanceFromLatLonInKm(52.295843, 4.9577432,
                                                   52.294629, 4.957973))
    distance_list.append(getDistanceFromLatLonInKm(40.7641718, -73.9567148,
                                                   40.764112, -73.955969))
    distance_list.append(getDistanceFromLatLonInKm(48.1508394, 11.5810441,
                                                   48.150839, 11.581044))
    distance_list.append(getDistanceFromLatLonInKm(44.6319678,10.943805,
                                                   44.688982, 10.663918))
    distance_list.append(getDistanceFromLatLonInKm(45.8238607,16.0077337,
                                                   51.524384, -0.073223))
    distance_list.append(getDistanceFromLatLonInKm(53.4630893,-2.2384236,
                                                   35.77601, 139.900834))
    distance_list.append(getDistanceFromLatLonInKm(43.0135332,-79.2667694,
                                                   43.013533, -79.266769))
    distance_list.append(getDistanceFromLatLonInKm(36.9914461,-122.0609145,
                                                   36.998541, -122.059815))
    distance_list.append(getDistanceFromLatLonInKm(-42.8778436,147.3287217,
                                                   -42.878589, 147.330118))
    distance_list.append(getDistanceFromLatLonInKm(13.1371586,77.497329,
                                                   13.134844, 77.4960078))

    num_correct = 0

    for d in distance_list:
        if d < 5:
            num_correct = num_correct + 1

    print('SUCCESS RATE: ' + str(num_correct / len(distance_list) * 100) + '%\n')
    print('MEAN DISTANCE: ' + str(statistics.mean(distance_list)) + ' km\n')
    print('1 STANDARD DEVIATION: ' + str(statistics.stdev(distance_list)) + ' km\n')

    print(str(distance_list))

if __name__ == "__main__":
    hitlist = fake_json('eFetch_hit.json')
    acclist = fake_json('eFetch_acc.json')

    accids = []

    sortedhits = sorted(hitlist, key=lambda k: k['MedlineCitation']['PMID'])

    for a in acclist:
        accids.append(a['MedlineCitation']['PMID'])

    for d in sortedhits:
        if d['MedlineCitation']['PMID'] in accids:
            print('\\textbf{' + d['MedlineCitation']['PMID'] + '}')
        else:
            print(d['MedlineCitation']['PMID'])
    # geocode.run(fake_json('eFetch_hit.json'))
    # doublecheck()
    # geocode_hit(fake_json('eFetch_hit.json'))
    # geocode_acc(fake_json('eFetch_acc.json'))
    # metamap.run(fake_json('mm_example.json'))
