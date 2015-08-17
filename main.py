import os
import json

import pymongo
from flask import Flask, render_template, request, jsonify
from pymongo import MongoClient
# from bson.objectid import ObjectId

import app.geocode, app.metamap, app.umls
import config
from app.metamap import MetaMap

# Flask setup
app = Flask(__name__, template_folder="./app/templates")

# MongoDB setup
client = MongoClient()
db = client.cached_results
pubmeddata = db.pubmeddata

# indexed by Pubmed ID in descending order
# db.pubmeddata.create_index([('MedlineCitation.PMID', pymongo.DESCENDING)])

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
    query = 'glioblastoma stem cells'
    places = geocode.get_data(query)
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    # app.run(debug=True)
    with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 'r') as datafile:
        results = json.load(datafile)
        print('JSON loaded')
        pubmeddata.insert(results[0])
        print('document inserted')
        cursor = db.pubmeddata.find_one( {'MedlineCitation.PMID': '00000000'} )
        if cursor:
            print('inserted and found woooo')
            print(cursor['MedlineCitation']['Article']['ArticleTitle'])
        else:
            print('booooo')
        datafile.close()

