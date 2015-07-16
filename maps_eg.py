from flask import Flask, render_template, request, geocode

app = Flask(__name__, template_folder=".")

def get_query():
    return request.form['text']    

@app.route("/")
def mapview():
    query = 'glioblastoma stem cells'
    # TODO: don't just pop an address!
    places = geocode.get_data(query)
    return render_template('maps_eg.html', places=places_dict)

if __name__ == "__main__":
    app.run(debug=True)

