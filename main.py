import os
import json
from multiprocessing import Process, Queue

import pymongo
from flask import Flask, render_template, request, jsonify
from pymongo import MongoClient

from app import metamap, geocode, umls

# Flask setup
app = Flask(__name__, template_folder="./app/templates")

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

def find_element(alist, pmid):
    for l in alist:
        if l['PMID'] is pmid:
            return l
    raise AttributeError('The expected PubMed ID is not in the list.')

# enter all results into the pubmeddata collection
def insert_into_db(results, concepts, places):
    for r in results:
        pmid = r['MedlineCitation']['PMID']
        r['concepts'] = find_element(concepts, pmid)
        r['placeids'] = find_element(places, pmid)
        insert_result = db.pubmeddata.insert_one(docs)

        # error logging
        if not insert_result.acknowledged:
            print('ERROR: not an acknowledged write operation.\n' +
                pmid + ' not inserted into db.')
        else: 
            print(insert_result.inserted_id + ' inserted')

@app.route("/")
def mapview():
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    # prevent circular references
    from app import citations

    query = 'glioblastoma stem cells'
    # docs = documents from MongoDB
    # results = new PubMed data
    (docs, results) = citations.start_search(query)
    
    # Multiprocessing using a queue
    # q = Queue()
    # m = Process(target=metamap.q_run, args=(results, q))
    # g1 = Process(target=geocode.q_run, args=(results, q))
    # g2 = Process(target=geocode.q_retrieve, args=(docs, q))
    
    # m.start()
    # g1.start()
    # g2.start()

    # block until item is available, so Process u can start
    # m_out = metamap.run(results)
    m_out = metamap.format_results(load_read_close('./tests/resources', 'metamap_output.txt').split('\n'))
    print(umls.run(m_out))
    # Retrieve UMLS hierarchy for both cached terms and terms from Metamap
    # combined = docs + m_out
    # print(umls.run(combined))
    
    # g1_out = q.get()
    # g2_out = q.get()
    
    # print(g1_out)
    # print(g2_out)

    # insert_into_db(results, m_out, gr_out)

    # app.run(debug=True)
