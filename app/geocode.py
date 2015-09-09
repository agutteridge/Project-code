import urllib
import urllib.error
import json
import os
import re

from app import config

# URLopen, read, decode, and turn into python object
def request(url):
    return (
        json.loads(
            urllib.request
            .urlopen(url)
            .read()
            .decode('UTF-8')
        )
    )

# Uses Google Places Web API to retrieve longitude and latitude corresponding to a place ID.
# Returns a python dict with the place name and coordinates
def get_latlong(placeid):
    url = ('https://maps.googleapis.com/maps/api/place/details/json?' +
            'placeid=' + placeid +
            '&key=' + config.maps_key +
            '&language=en')

    result = request(url)

    # retain structure for continuity between data from cache and API
    return {
        'name': result['result']['name'],
        'geometry': {
            'location': result['result']['geometry']['location']
        }
    }

# Uses Google Places Web API to match the address string to a real-world location.
# Returns a list 
def get_location(num, address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    # establishment is Google default
    url = ('https://maps.googleapis.com/maps/api/place/textsearch/json?' + 
        encoded_query +
        '&key=' + config.maps_key +
        '&types=university|hospital|establishment' +
        '&language=en') 

    place_options = []

    try:
        place_options = request(url)['results']
    except (urllib.error.HTTPError,
            urllib.error.URLError) as e:
        print(e.reason)
        print(url)

    result = dict()

    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        datafile.write('Formatted address (format' + str(num) + '): ' + address + '\n')
        if place_options:
            datafile.write('Output address: ' +  place_options[0]['name'] + '\n')
            datafile.write('Coordinates: ' + str(place_options[0]['geometry']['location']['lat']) + 
                ", " + str(place_options[0]['geometry']['location']['lng']) + '\n')
            result = place_options[0] # only return first result
        elif num is 3:
            datafile.write('Output address: NONE\n')
        datafile.close()

    return result

# Email addresses are removed to improve success of geocoding using Google Places API
def remove_email(address):
    # regex for removing nonwhitespace@[alphanum-.]+
    without_chars = re.sub('[\S]+[@][\w.-]+', '', address) 
    pattern = re.compile('(email|address|electronic)', re.IGNORECASE)
    without_email = pattern.sub('', without_chars)
    result = without_email.split(',')
    return result

# Removes lines from address that refer to a department
def remove_dept(address):
    address_lines = address.split(',')
    results = []
    for a in address_lines:
        temp_a = a.upper()
        if 'DEPT' not in temp_a and 'DEPARTMENT' not in temp_a:
            results.append(a.strip()) # remove trailing whitespace
    return results

# Remove email only
def format0(address):
    without_email = remove_email(address)
    formatted_address = ','.join(without_email).strip()
    return formatted_address

# Remove first line
def format1(address):
    without_email = remove_email(address)
    if len(without_email) > 1:
        return ','.join(without_email[1:len(without_email)]).strip()
    else:
        return without_email.strip()

# Second and third lines only
def format2(address):
    without_email = remove_email(address)
    if len(without_email) > 2:
        return ', '.join(without_email[1:3]).strip()
    else: # Less than 3 lines, use last line only
        return without_email[-1:].strip()

def format3(address):
    without_email = remove_email(address)
    if len(without_email) > 2:
        return ', '.join([without_email[1], without_email[-1]]).strip()
    else: # Less than 3 lines, use first line only
        return without_email[0].strip()

# Returns a list of dicts, each with a different address.
# The list of strings with alphanumeric chars prevents duplication of
# addresses with minor punctuation differences
def unique_addresses(addresses):
    unique_list = []

    for a in addresses:
        address = a['address']
        # exclude all non-alphanumeric chars, upper case
        alphanumeric = re.sub('[\W]', '', address).upper()
        added = False

        for u in unique_list:
            if u['alphanumeric'] == alphanumeric:
                u['i_nums'].append(a['i_num'])

                if a['pmid'] not in u['pmids']:
                    u['pmids'].append(a['pmid'])
            
                added = True

        if not added:
            unique_list.append({
                'alphanumeric': alphanumeric,
                'address': address,
                'pmids': [a['pmid']],
                'i_nums': [a['i_num']]
            })

    return unique_list

# Returns a list of all addresses, sorted into dicts with the corresponding PubMed ID
def list_all_addresses(papers):
    all_addresses = []

    i_num = 0 # assign unique ID to each paper
    for paper in papers:
        pmid = str(paper['MedlineCitation']['PMID'])

        for author in paper['MedlineCitation']['Article']['AuthorList']:
            for place in author['AffiliationInfo']:
                individual_addresses = place['Affiliation'].split(';')

                for i in individual_addresses:
                    formatted_address = format0(i) # Remove email from all addresses
                    all_addresses.append({
                        'pmid': pmid, 
                        'address': formatted_address,
                        'i_num': i_num })
                    i_num = i_num + 1

    return all_addresses

def remove_address(i_num, all_addresses):
    remaining_addresses = []
    
    for a in all_addresses:
        if a['i_num'] is not i_num:
            remaining_addresses.append(a)

    return remaining_addresses

def start_geocoding(all_addresses):
    geocoded_addresses = []
    unique0 = unique_addresses(all_addresses)

    # Round 0, email removed
    for u in unique0:
        result0 = get_location(0, u['address'])

        if result0:
            for p in u['pmids']:
                geocoded_addresses.append({'PMID': p, 'place': result0 })
            for i in u['i_nums']:
                all_addresses = remove_address(i, all_addresses) # mutate original list of addresses

    # Round 1, first line removed
    formatted1 = []

    for a in all_addresses:
        formatted_address = format1(a['address'])
        formatted1.append({'pmid': a['pmid'], 
            'address': formatted_address,
            'i_num': a['i_num'] })

    unique1 = unique_addresses(formatted1)

    for u in unique1:
        result1 = get_location(1, u['address'])

        if result1:
            for p in u['pmids']:
                geocoded_addresses.append({'PMID': p, 'place': result1 })
            for i in u['i_nums']:
                all_addresses = remove_address(i, all_addresses)

    # Round 2, second and third lines only
    formatted2 = []

    for a in all_addresses:
        formatted_address = format2(a['address'])
        formatted2.append({'pmid': a['pmid'], 
            'address': formatted_address,
            'i_num': a['i_num'] })

    unique2 = unique_addresses(formatted2)

    for u in unique2:
        result2 = get_location(2, u['address'])

        if result2:
            for p in u['pmids']:
                geocoded_addresses.append({'PMID': p, 'place': result2 })
            for i in u['i_nums']:
                all_addresses = remove_address(i, all_addresses)

    # Round 2, second and third lines only
    formatted3 = []

    for a in all_addresses:
        formatted_address = format3(a['address'])
        formatted3.append({'pmid': a['pmid'], 
            'address': formatted_address,
            'i_num': a['i_num'] })

    unique3 = unique_addresses(formatted3)

    for u in unique3:
        result3 = get_location(3, u['address'])

        if result3:
            for p in u['pmids']:
                geocoded_addresses.append({'PMID': p, 'place': result3 })
            for i in u['i_nums']:
                all_addresses = remove_address(i, all_addresses)

    print(str(len(all_addresses)) + " addresses not geocoded.")
    for leftover in all_addresses:
        geocoded_addresses.append({'PMID': leftover['pmid'], 'place': dict()})

    return geocoded_addresses

def run(results):
    print('in geocode.run')

    # dicts with PMID and original addresses
    all_addresses = list_all_addresses(results)
    results = start_geocoding(all_addresses)
    for_cache = []

    # Sorting by PMID for cache
    for r in results:
        if r['place']:
            added = False
            
            for c in for_cache:
                if c['PMID'] == r['PMID']:
                    c['placeids'].append(r['place']['place_id'])
                    added = True

            if not added:
                for_cache.append({'PMID': r['PMID'],
                    'placeids': [r['place']['place_id']]})

    print('returning from geocode.run')
    return {
        'results': results,
        'for_cache': for_cache
    }

    # Returns a list of dicts
def retrieve(docs):
    print('in geocode.retrieve')
    results = []

    for d in docs:
        pmid = d['MedlineCitation']['PMID']   
        
        for p in d['placeids']:
            results.append( {'PMID': pmid, 'place': get_latlong(p)} )

    print('returning from geocode.retrieve')        
    return results
