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

def insert_into_db(docs):
    # enter all results into the pubmeddata collection
    insertmanyresult = db.pubmeddata.insert_many(docs)

    # error logging
    if not insertmanyresult.acknowledged:
        print('ERROR: not an acknowledged write operation.')

    if len(insertmanyresult.inserted_ids) is not len(docs):
        print('ERROR: The following documents were not inserted:')
        for d in docs:
            if d not in insertmanyresult.inserted_ids:
                print(d)
        print('not inserted.')

@app.route("/")
def mapview():
    return render_template('maps_eg.html', places=places)

if __name__ == "__main__":
    # prevent circular references
    from app.citations import start_search
    from app import metamap

    query = 'glioblastoma stem cells'
    # docs = documents from MongoDB
    # results = PubMed data
    (docs, results) = start_search(query)
    
    parent_conn_m, child_conn_m = Pipe()
    # Multiprocessing
    m = Process(target=metamap.pipe_in, args=(results, child_conn_m))
    # g1 = Process(target=geocode.run, args=results)
    # g2 = Process(target=geocode.get_placeids, args=docs)
    
    m.start()
    # g1.start()
    # g2.start()

    m_out = parent_conn_m.recv()
    m.join()
    # g1_out = parent_conn.recv()
    # g1.join()
    # g2_out = parent_conn.recv()
    # g2.join()
    print(m_out)

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

