import os
import json

import pymongo
from flask import Flask, render_template, request, jsonify
from pymongo import MongoClient


# from app.metamap import MetaMap

# Flask setup
app = Flask(__name__, template_folder="./app/templates")

# MongoDB setup
client = MongoClient()
db = client.cached_results

# indexed by Pubmed ID in descending order
db.pubmeddata.create_index([('MedlineCitation.PMID', pymongo.DESCENDING)])

# indexed by place ID in ascending order
# db.placeids.create_index(('placeID', pymongo.ASCENDING))

# module level functions
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
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    # prevent circular references
    from app.citations import start_search
    query = 'glioblastoma stem cells'
    places = start_search(query)
    # need to also return list of IDs or some kind of indicator for papers that are just now being added to the cache,
    # so that metamap CUIs and place IDs can also be added to the DB.
    print(places)
    # app.run(debug=True)
