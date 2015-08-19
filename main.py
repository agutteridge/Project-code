import os
import json
from multiprocessing import Process, Pipe

import pymongo
from flask import Flask, render_template, request, jsonify
from pymongo import MongoClient

# Flask setup
app = Flask(__name__, template_folder="./app/templates")

# MongoDB setup
client = MongoClient()
db = client.cached_results

# if not db.pubmeddata.list_indexes():
    # indexed by Pubmed ID in descending order
    # db.pubmeddata.create_index([('MedlineCitation.PMID', pymongo.DESCENDING)])

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
        r['places'] = find_element(places, pmid)[]
        insert_result = db.pubmeddata.insert_one(docs)

        # error logging
        if not insert_result.acknowledged:
            print('ERROR: not an acknowledged write operation.\n' +
                pmid + ' not inserted into db.')

@app.route("/")
def mapview():
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    # prevent circular references
    from app import citations
    from app import metamap

    query = 'glioblastoma stem cells'
    # docs = documents from MongoDB
    # results = new PubMed data
    (docs, results) = citations.start_search(query)
    
    # Multiprocessing using a queue
    q = Queue()
    m = Process(target=metamap.q_run, args=(results, q))
    gr = Process(target=geocode.q_run, args=(results, q))
    # g2 = Process(target=geocode.get_placeids, args=docs)
    
    m.start()
    gr.start()
    # g2.start()

    m_out = m.get()
    gr_out = fr.get()
    m.join()
    # g1_out = parent_conn.recv()
    # g1.join()
    # g2_out = parent_conn.recv()
    # g2.join()
    print(m_out)

    insert_into_db(results, m_out, gr_out)

    # app.run(debug=True)

# # testing..
#     # with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 'r') as datafile:
#     #     mm = MetaMap()
#     #     mm.run(json.load(datafile))
#     #     datafile.close()

#     test_txt = load_read_close('./tests/resources', 'metamap_output.txt').split('\n')
#     paper_terms = app.metamap.format_results(test_txt)
#     results = app.umls.run_with_data(paper_terms)
#     for r in results:
#     	print(r)

