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
def get_location(address):
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

    #logging
    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        datafile.write('Input address: ' + address + '\n')
        if place_options:
            datafile.write('Output address: ' +  place_options[0]['name'] + '\n')
            result = place_options[0] # only return first result
        else:
            datafile.write('Output address: NONE\n')
        datafile.close()

    return result

# Email addresses are removed to improve success of geocoding using Google Places API
def remove_email(address):
    # regex for removing nonwhitespace@[alphanum-.]+
    result = re.sub('[\S]+[@][\w.-]+', '', address) 
    reg = re.compile('(email|address|electronic)', re.IGNORECASE)
    result = reg.sub('', result)
    # TODO: make sure comma is also deleted otherwise last thing will be just crap
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

# Removes first address lines for addresses over ???????? parts long (comma separated)
def format_address(address):
    without_email = remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 1:
        return ','.join(address_lines[1:len(address_lines)])
    else:
        return ''

# Returns a set of strings, each with a different address.
# The set of strings with alphanumeric chars prevents duplication of
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
                formatted_address = format_address(f)
                # exclude all non-alphanumeric chars, upper case
                alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                if alphanumeric != '' and alphanumeric not in alphanumeric_addresses and formatted_address:
                    result.append(formatted_address)
                    alphanumeric_addresses.add(alphanumeric)
    return result

def run(results):
    print('in geocode.run')
    result_list = []
    for_cache = []

    for paper in results:
        pmid = str(paper['MedlineCitation']['PMID'])
        author_list = paper['MedlineCitation']['Article']['AuthorList']
        addresses = (unique_addresses(author_list))
        
        placeids = []

        for address in addresses:
            place = get_location(address)
            if place:
                result_list.append( {'PMID': pmid, 'place': place } )
                placeids.append(place['place_id'])

        for_cache.append({'PMID': pmid,
            'placeids': placeids})

    print('returning from geocode.run')
    return {
        'results': result_list,
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
