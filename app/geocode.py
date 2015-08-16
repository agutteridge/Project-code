import urllib
import sys
import json
import os

import config
from app import citations

def get_location(address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    full_url = ('https://maps.googleapis.com/maps/api/place/textsearch/json?' +
                encoded_query +
                '&key=' + config.maps_key +
                '&types=' + 'university|hospital|establishment') # establishment is Google default

    result = urllib.request.urlopen(full_url)
    return result.read()

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

# write a function for returning all useful data, parsed nicely as JSON
# relevant fields:
# all authors' names
# Institutes of all authors (names only)
# Institutes of first & last author (both long and lat and names)
# journal
# paper title
# year
# MeSH terms
# keywords
