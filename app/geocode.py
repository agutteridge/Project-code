import urllib
import sys
import json
import app.config
import app.citations

def get_location(address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    full_url = 'https://maps.googleapis.com/maps/api/place/textsearch/json?' + encoded_query + '&key=' + config.maps_key
    result = urllib.request.urlopen(full_url)
    return result.read()

def get_data(query):
    addresses = citations.get_addresses(query)

    successful_address = 'none'
    for i in range(0, len(addresses)):
        if addresses[i]:
            successful_address = addresses[i]
            break

    print("address to be geocoded: " + successful_address)
    places_bytes = get_location(successful_address)
    places_str = places_bytes.decode('UTF-8')
    # parse into long and lat only
    return places_str

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
