from flask import Flask, render_template
from flask_googlemaps import GoogleMaps
from flask_googlemaps import Map
import urllib
import sys
import json

app = Flask(__name__, template_folder=".")
GoogleMaps(app)
api_key = sys.argv[1]

def get_location(address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    full_url = 'https://maps.googleapis.com/maps/api/place/textsearch/json?' + encoded_query + '&key=' + api_key
    result = urllib.request.urlopen(full_url)
    return result.read()

@app.route("/")
def mapview():
    # getting location
    results = get_location('Birkbeck, University of London')
    places_str = results.decode('UTF-8')
    print(places_str)
    places_dict = json.loads(places_str)

    lat = places_dict['results'][0]['geometry']['location']['lat']
    print(lat)
    lng = places_dict['results'][0]['geometry']['location']['lng']
    print(lng)

    # creating a map in the view
    mymap = Map(
        identifier="view-side",
        lat = lat,
        lng = lng,
        markers=[(lat, lng)]
    )
    return render_template('maps_eg.html', mymap=mymap)

if __name__ == "__main__":
    app.run(debug=True)
