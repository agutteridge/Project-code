import os
import datetime
import pymongo
from pymongo import MongoClient

# MongoDB setup
client = MongoClient()
db = client.cached_results
pubmeddata = db['pubmeddata']

def find_element(placeids_or_concepts, pmid, string):
    for pc in placeids_or_concepts:
        if pc['PMID'] == pmid:
            return pc[string]

    print(str(pmid) + ' has no ' + string)
    return []

# enter all results into the pubmeddata collection
def insert_into_db(results, concepts, places):
    all_docs = []

    for r in results:
        pmid = str(r['MedlineCitation']['PMID'])
        r['concepts'] = find_element(concepts, pmid, 'concepts')
        r['placeids'] = find_element(places, pmid, 'placeids')

        #logging
        with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
            datafile.write('\n' + pmid + ':\n\t' + str(len(r['concepts'])) + ' concepts.\n') 
            
            for c in r['concepts']:
                datafile.write('\t\t' + str(c) + '\n')

            datafile.write('\t' + str(len(r['placeids'])) + ' place IDs.\n')
            
            for p in r['placeids']:
                datafile.write('\t\t' + str(p) + '\n')

            datafile.write("----------------------------------------------------------------------------\n")

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

def retrieve(id_list):
    new_ids = []
    cached_docs = []

    for i in id_list:
        # search for PMID in pubmeddata mongoDB collection
        cursor = db.pubmeddata.find_one( {'MedlineCitation.PMID': i} )

        if not cursor: # if paper not already cached
            new_ids.append(i)
        else:
            cached_docs.append(cursor)

    return {
        'new': new_ids,
        'old': cached_docs
    }
