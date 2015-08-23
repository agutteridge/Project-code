import os
import json
from multiprocessing import Process, Queue

import pymongo
from flask import Flask, render_template, request, jsonify
from pymongo import MongoClient

from app import metamap, geocode, umls

# Flask setup
app = Flask(__name__, template_folder='./app/templates')

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
        print(r)
        all_docs.append(r)

    # batch all insertions
    insert_result = db.pubmeddata.insert_many(all_docs)

    if not insert_result.acknowledged:
        print('ERROR: not an acknowledged write operation.')
    else: 
        print(str(len(insert_result.inserted_ids)) + ' inserted')

@app.route("/")
def mapview():
    # return render_template('maps_eg.html', places=places)
    return render_template('d3_eg.html')

@app.route("/data")
def return_data():
    with open(os.path.join('./app/static', 'umls_output.json'), 'r') as datafile:
        obj = json.load(datafile)
        datafile.close()
        return json.dumps(obj)

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
    
    (front, back) = geocode.run(results)

    # m.start()
    # g1.start()
    # g2.start()

    # block until item is available, so Process u can start
    m_out = metamap.run(results)
    # m_out = metamap.format_results(load_read_close('./tests/resources', 'metamap_output.txt').split('\n'))
    # umls.run(m_out))

    # Retrieve UMLS hierarchy for both cached terms and terms from Metamap
    combined = docs + m_out
    umls.run(combined)
    
    # g1_out = q.get()
    # g2_out = q.get()
    
    # print(g1_out)
    # print(g2_out)

    if len(m_out) > 0 or len(back) > 0:
        insert_into_db(results, m_out, back)
    else:
        print('all cached baby')

    app.run(debug=True)
