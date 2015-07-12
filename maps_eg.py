from flask import Flask, render_template, request
import urllib, sys, json, printing_affiliations, config

app = Flask(__name__, template_folder=".")

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
    # TODO: don't just pop an address!
    random_address = printing_affiliations.get_addresses(query).pop()
    print(random_address)
    places_bytes = get_location(random_address)
    places_str = places_bytes.decode('UTF-8')
    places_dict = json.loads(places_str)

    if (places_dict['results']): 
        print(places_dict['results'])
        return render_template('maps_eg.html', places=places_dict)
    else:
        # london as default
        return render_template('maps_eg.html', 
            places={"results" : [{"geometry" : {"location" : {"lat" : 51.508742, "lng" : -0.120850}}}]})

if __name__ == "__main__":
    app.run(debug=True)

