from flask import Flask, render_template, request, jsonify
import geocode

app = Flask(__name__, template_folder=".")

def get_query():
    return request.form['text']    

@app.route("/")
def mapview():
    query = 'glioblastoma stem cells'
    places = geocode.get_data(query)
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    app.run(debug=True)

