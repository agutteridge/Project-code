import concurrent.futures

import os
import datetime
import json

from app import citations, geocode, metamap, umls, cache

def search_for(query):
    #logging
    with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
        datafile.write('BEGIN ' + str(datetime.datetime.today()) + '\n')
        datafile.write('Search term: ' + query + '\n')
        datafile.close()

    citations_results = citations.start_search(query)

    if citations_results:
        #logging
        with open(os.path.join('./app/static', 'log.txt'), 'a') as datafile:
            datafile.write(str(len(citations_results['docs'])) + ' documents retrieved from database.\n')
            datafile.write(str(len(citations_results['results'])) + ' documents fetched from PubMed.\n')
            datafile.close()

        executor = concurrent.futures.ThreadPoolExecutor(max_workers=3)
        m = executor.submit(metamap.run, citations_results['results'])
        g = executor.submit(geocode.run, citations_results['results'])
        r = executor.submit(geocode.retrieve, citations_results['docs'])

        metamap_results = m.result()
        combined = citations_results['docs'] + metamap_results

        u = executor.submit(umls.run, combined)

        geocode_run_results = g.result() 
        places = geocode_run_results['results'] + r.result()
        concepts = u.result()

        if len(metamap_results) > 0 or len(geocode_run_results['for_cache']) > 0:
            cache.insert_into_db(citations_results['results'], metamap_results, geocode_run_results['for_cache'])
        else:
            print('all docs retrieved from cache')

        return json.dumps({ 'concepts': concepts, 'places': places, 'papers': citations_results['formatted'] })
    else:
        return json.dumps('No results found')