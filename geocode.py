import urllib, sys, json, config, citations

def get_location(address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    full_url = 'https://maps.googleapis.com/maps/api/place/textsearch/json?' + encoded_query + '&key=' + config.maps_key
    result = urllib.request.urlopen(full_url)
    return result.read()

def get_data(query):
	# instead of pop..?
	random_address = citations.get_addresses(query).pop()
    places_bytes = get_location(random_address)
    places_str = places_bytes.decode('UTF-8')
    places_dict = json.loads(places_str)

    # parse into long and lat only
    return places_dict

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
