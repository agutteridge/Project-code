import os
import json
from multiprocessing import Process, Queue
import datetime
import time

import pymongo
from flask import Flask, render_template, request, jsonify
from pymongo import MongoClient

from app import metamap, geocode, umls

# Flask setup
app = Flask(__name__, template_folder='./app/templates',
                      static_folder='./app/static')

# MongoDB setup
client = MongoClient()
db = client.cached_results
pubmeddata = db['pubmeddata']

def load_read_close(path, filename):
    with open(os.path.join(path, filename), 'r') as datafile:
        txt = datafile.read()
        datafile.close()
        return txt

# TODO
def get_query():
    return request.form['text']    

def find_element(results, pmid, placeids_or_concepts):
    for l in results:
        if l['MedlineCitation']['PMID'] == pmid:
            return l[placeids_or_concepts]
    raise AttributeError('The expected PubMed ID is not in the list.')

# enter all results into the pubmeddata collection
def insert_into_db(results, concepts, places):
    all_docs = []

    for r in results:
        pmid = r['MedlineCitation']['PMID']
        r['concepts'] = find_element(concepts, pmid, 'concepts')
        r['placeids'] = find_element(places, pmid, 'placeids')

        #logging
        with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
            datafile.write(pmid + ': ' + str(len(r['concepts'])) + 
                ' concepts, ' + str(len(r['placeids'])) + ' place IDs.\n')
            datafile.close()

        all_docs.append(r)

    # batch all insertions
    insert_result = db.pubmeddata.insert_many(all_docs)

    if not insert_result.acknowledged:
        print('ERROR: not an acknowledged write operation.')
    else: 
        with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
            datafile.write(str(len(insert_result.inserted_ids)) + ' inserted:\n')
            for i in insert_result.inserted_ids:
                datafile.write(str(i) + '\n')
            datafile.close()

@app.route("/")
def index():
    # get another route through which JSON can be sent woooooooooo (maybe same structure for both D3 and GMaps??)
    return render_template('index.html')

@app.route("/data")
def return_data():
    search_term = request.args.get('search_term', '')
    return search_for(search_term)

def search_for(query):
    #logging
    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        now = datetime.datetime.today()
        datafile.write('BEGIN ' + str(now) + '\n')
        datafile.write('Search term: ' + query + '\n')
        datafile.close()

    # prevent circular references
    from app import citations

    (docs, results) = citations.start_search(query)
    
    #logging
    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        datafile.write(str(len(docs)) + ' documents retrieved from database.\n')
        datafile.write(str(len(results)) + ' documents fetched from PubMed.\n')
        datafile.close()

    m_out = metamap.run(results)
    combined = docs + m_out
    concepts = umls.run(combined)

    (front1, back) = geocode.run(results)
    front2 = geocode.retrieve(docs)
    places = front1 + front2

    if len(m_out) > 0 or len(back) > 0:
        insert_into_db(results, m_out, back)
    else:
        print('all docs retrieved from cache')

    return json.dumps({ 'concepts' : concepts, 'places' : places })


if __name__ == "__main__":
    app.run(debug=True)
