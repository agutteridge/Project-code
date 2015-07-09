from flask import Flask, render_template, request
import urllib, sys, json, printing_affiliations, config

app = Flask(__name__, template_folder=".")
# api_key = sys.argv[1]

def get_location(address):
    query_string = {'query': address}
    encoded_query = urllib.parse.urlencode(query_string)
    full_url = 'https://maps.googleapis.com/maps/api/place/textsearch/json?' + encoded_query + '&key=' + config.maps_key
    result = urllib.request.urlopen(full_url)
    return result.read()

def get_query():
    return request.form['text']    

@app.route("/")
def mapview():
    query = 'glioblastoma stem cells'
    random_address = printing_affiliations.get_addresses(query).pop()
    # print(random_address)
    # # getting location
    # places_bytes = get_location(random_address)
    # places_str = places_bytes.decode('UTF-8')
    # # print(places_str)
    # places_dict = json.loads(places_str)

    # return render_template('maps_eg.html', places=places_dict)
    return random_address

if __name__ == "__main__":
    app.run(debug=True)

