import urllib
import sys
import json
import os

import config
from app import citations

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
                formatted_address = _format_address(f)
                # exclude all non-alphanumeric chars, upper case
                alphanumeric = re.sub('[\W]', '', formatted_address).upper()
                if alphanumeric != '' and alphanumeric not in alphanumeric_addresses:
                    result.append(formatted_address)
                    print("New address: " + formatted_address)
                    alphanumeric_addresses.add(alphanumeric)
    print("List of unique addresses: " + str(result))
    return result

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
    without_email = _remove_email(address)
    address_lines = without_email.split(',')
    if len(address_lines) > 2:
        address_lines = address_lines[-2:] # last 2 elements of the list
    result = (','.join(address_lines))
    return result


# Uses Google Places Web API to match the address string to a real-world location.
def get_location(address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    full_url = ('https://maps.googleapis.com/maps/api/place/textsearch/json?' +
                encoded_query +
                '&key=' + config.maps_key +
                '&types=' + 'university|hospital|establishment') # establishment is Google default

    result = urllib.request.urlopen(full_url)
    return result.read()

# rename to run or something ???
def get_addresses(papers):
    result_list = []

    for paper in papers:
        author_list = paper['MedlineCitation']['Article']['AuthorList']
        result_list = result_list + (unique_addresses(author_list))

    return result_list

# can be replaced with tests/resources/geocode_output.json
def get_data(query):
    addresses = citations.get_addresses(query)

    points = []
    for a in addresses:
        places_bytes = get_location(a)
        places_str = places_bytes.decode('UTF-8')
        points.append(json.loads(places_str))
    # parse into long and lat only

    return json.dumps(points)