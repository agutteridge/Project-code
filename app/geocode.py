import urllib
import sys
import json
import os
import re

import config

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
            '&key=' + config.maps_key)

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
    url = ('https://maps.googleapis.com/maps/api/place/textsearch/json?' + 
        encoded_query +
        '&key=' + config.maps_key + 
        '&types=university|hospital|establishment') # establishment is Google default

    place_options = request(url)['results']

    if place_options:
        return place_options[0] # only return first result
    else: 
        return dict()

# Email addresses are removed to improve success of geocoding using Google Places API
def remove_email(address):
    # regex for removing nonwhitespace@[alphanum-.]+
    result = re.sub('[\S]+[@][\w.-]+', '', address) 
    reg = re.compile('(email|address|electronic)', re.IGNORECASE)
    result = reg.sub('', result)
    # TODO: make sure comma is also deleted otherwise last thing will be just crap
    return result

# Removes first address lines for addresses over ???????? parts long (comma separated)
def format_address(address):
    without_email = remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        address_lines = address_lines[-2:] # last 2 elements of the list
    result = (','.join(address_lines))
    return result

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
                if alphanumeric != '' and alphanumeric not in alphanumeric_addresses:
                    result.append(formatted_address)
                    alphanumeric_addresses.add(alphanumeric)
    return result

def q_run(results, q):
    print('q_run start in geocode')
    q.put(run(results))
    print('q_run in geocode done')

def q_retrieve(docs, q):
    print('q_retrieve start in geocode')
    q.put(retrieve(docs))
    print('q_retrieve in geocode done')

def run(results):
    result_list = []
    all_cache = []

    for paper in results:
        pmid = str(paper['MedlineCitation']['PMID'])
        author_list = paper['MedlineCitation']['Article']['AuthorList']
        addresses = (unique_addresses(author_list))
        address_list = []
        
        for address in addresses:
            address_list.append(get_location(address))

        result_list.append({'PMID': pmid, 'places': address_list})

        placeids = []
        
        for a in address_list:
            if a:
                placeids.append(a['place_id'])

        all_cache.append({'MedlineCitation': {'PMID': pmid},
            'placeids': placeids})

    return (result_list, all_cache)

    # Returns a list of dicts
def retrieve(docs):
    results = []

    for d in docs:
        doc_places = dict()
        doc_places['PMID'] = d['MedlineCitation']['PMID']   
        placeid_list = []
        
        for p in d['placeids']:
            placeid_list.append(get_latlong(p))
        
        doc_places['results'] = placeid_list
        results.append(doc_places)

    return results
