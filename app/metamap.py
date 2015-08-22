# The Java .jar that calls the Medical Text Indexer is invoked from this file

import subprocess
import datetime
import os
import time
import re

import config

# Creates a batch ID for .txt files to be used as input for MetaMap
def create_batch_id():
    # removes last 4 digits so ms is 2 s.f.
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-4] 

# Formats MetaMap results into a list of dicts, each with a PubMed ID and a list of 
# lists of concepts for that paper
def format_results(results):
    all_concepts = []
    for r in results:
        concept = r.split('|')
        PMID = concept[0]
        PMID_non_digits = re.sub('[\d]', '', PMID)
        rest = concept[1:]

        # filtering out empty strings and those with non-digits
        if len(PMID_non_digits) is 0 and len(PMID) is not 0:
            if not all_concepts or all_concepts[0]['PMID'] != PMID:
                # keeping structure
                all_concepts = [
                    {
                    'MedlineCitation': 
                        {
                        'PMID': PMID
                        }, 
                    'concepts': []
                    }
                ] + all_concepts # prepend
            all_concepts[0]['concepts'].append(rest) # adds concept to list of concepts for that paper
    print(all_concepts)
    return all_concepts

# creates a text file for MetaMap to use as a source
def write_file(filename, results):
    batch = open(os.path.join('app/static', filename), 'w') # default: unbuffered
    for i in range(0, len(results)):
        # exclude ASCII characters to avoid MetaMap errors
        ASCII_title = results[i]['MedlineCitation']['Article']['ArticleTitle'].encode('ascii', 
            errors='ignore').decode('UTF-8')
        ASCII_abstract = results[i]['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode('ascii', 
            errors='ignore').decode('UTF-8')
        batch.write('UI  - ' + results[i]['MedlineCitation']['PMID'] + '\n' + 
            'TI  - ' + ASCII_title + '\n' +
            'AB  - ' + ASCII_abstract + 
            '\n\n')
    batch.close()
    return True

def q_run(results, q):
  print('q_run start in metamap')
  q.put(run(results))
  print('q_run in metamap done')

def run(results):
    batch_id = create_batch_id()
    filename = batch_id + '.txt'
    write_file(filename, results)

    mp = config.metamap_path

    # shell args for running MetaMapCaller.class
    popen_args = ['java',
                  '-cp',
                  mp + '/classes:' +
                  mp + '/lib/skrAPI.jar:' +
                  mp + '/lib/commons-logging-1.1.1.jar:' +
                  mp + '/lib/httpclient-cache-4.1.1.jar:' +
                  mp + '/lib/httpcore-nio-4.1.jar:' +
                  mp + '/lib/httpclient-4.1.1.jar:' +
                  mp + '/lib/httpcore-4.1.jar:' +
                  mp + '/lib/httpmime-4.1.1.jar:.',
                  'MetaMapCaller',
                  '../../static/' + filename, 
                  config.un, # username for MetaMap
                  config.pwd, # password for MetaMap
                  config.email, # email for MetaMap
                  '-y'] 

    p = subprocess.Popen(popen_args, cwd='app/java/bin/', stdout=subprocess.PIPE)

    terms_list = []

    while 1:
        term = p.stdout.readline()
        if not term and p.returncode is not None:
            break
        terms_list.append(term.decode('UTF-8'))
        p.poll()
    print("done %d" % p.returncode)

    return format_results(terms_list)
