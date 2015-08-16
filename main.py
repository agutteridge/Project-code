from flask import Flask, render_template, request, jsonify
from flask.ext.mongoengine import MongoEngine
import os
import json

from app import geocode
import app.metamap
import app.umls
from app.metamap import MetaMap
import config

app = Flask(__name__, template_folder="./app/templates")

app.config['MONGODB_SETTINGS'] = {'DB': 'cached_results'}
app.config['SECRET_KEY'] = config.secret

db = MongoEngine(app)

def load_read_close(path, filename):
    with open(os.path.join(path, filename), 'r') as datafile:
        txt = datafile.read()
        datafile.close()
        return txt

# # testing..
# if __name__ == "__main__":
#     # with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 'r') as datafile:
#     #     mm = MetaMap()
#     #     mm.run(json.load(datafile))
#     #     datafile.close()

#     test_txt = load_read_close('./tests/resources', 'metamap_output.txt').split('\n')
#     paper_terms = app.metamap.format_results(test_txt)
#     results = app.umls.run_with_data(paper_terms)
#     for r in results:
#     	print(r)


def get_query():
    return request.form['text']    

@app.route("/")
def mapview():
    query = 'glioblastoma stem cells'
    places = geocode.get_data(query)
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    app.run(debug=True)

