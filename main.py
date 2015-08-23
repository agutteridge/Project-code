import os
import json
from multiprocessing import Process, Queue
import datetime

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

json_for_d3 = json.dumps(dict())

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
    return render_template('d3_eg.html')

@app.route("/data")
def return_data():
    search_term = request.args.get('search_term', '')
    return search_for(search_term)

def search_for(query):
    #logging
    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        time = datetime.datetime.today()
        datafile.write('BEGIN ' + str(time) + '\n')
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

    # Multiprocessing using a queue
    # q = Queue()
    # m = Process(target=metamap.q_run, args=(results, q))
    # g1 = Process(target=geocode.q_run, args=(results, q))
    # g2 = Process(target=geocode.q_retrieve, args=(docs, q))
    
    (front, back) = geocode.run(results)
    g_out = geocode.retrieve(docs)

    # m.start()
    # g1.start()
    # g2.start()

    # block until item is available, so Process u can start
    m_out = metamap.run(results)
    # m_out = metamap.format_results(load_read_close('./tests/resources', 'metamap_output.txt').split('\n'))
    # umls.run(m_out))

    # Retrieve UMLS hierarchy for both cached terms and terms from Metamap
    combined = docs + m_out

    # g1_out = q.get()
    # g2_out = q.get()
    
    # print(g1_out)
    # print(g2_out)

    if len(m_out) > 0 or len(back) > 0:
        insert_into_db(results, m_out, back)
    else:
        print('all cached baby')

    return umls.run(combined)


if __name__ == "__main__":
    app.run(debug=True)
