import os
import datetime
import json

from app import citations, geocode, metamap, umls, cache

def search_for(query):
    start = ''

    #logging
    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        now = datetime.datetime.today()
        start = str(now)
        datafile.write('BEGIN ' + start + '\n')
        datafile.write('Search term: ' + query + '\n')
        datafile.close()

    citations_results = citations.start_search(query)

    if citations_results:
        #logging
        with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
            datafile.write(str(len(citations_results['docs'])) + ' documents retrieved from database.\n')
            datafile.write(str(len(citations_results['results'])) + ' documents fetched from PubMed.\n')
            datafile.close()

        metamap_results = metamap.run(citations_results['results'])
        combined = citations_results['docs'] + metamap_results
        concepts = umls.run(combined)

        geocode_run_results = geocode.run(citations_results['results'])
        geocode_retrieve_results = geocode.retrieve(citations_results['docs'])
        places = geocode_run_results['results'] + geocode_retrieve_results

        if len(metamap_results) > 0 or len(geocode_run_results['for_cache']) > 0:
            cache.insert_into_db(citations_results['results'], metamap_results, geocode_run_results['for_cache'])
        else:
            print('all docs retrieved from cache')

        end(start)

        return json.dumps({ 'concepts': concepts, 'places': places, 'papers': citations_results['formatted'] })
    else:
        end(start)
        return 'No items found.'

def end(start):
    #logging
    with open(os.path.join('./app/static', 'speedtest.txt'), 'a') as datafile:
        now = datetime.datetime.today()
        datafile.write('BEGIN ' + start + '\n')
        datafile.write('END ' + str(now) + '\n')
        datafile.close()    