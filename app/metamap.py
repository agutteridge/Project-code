# The Java .jar that calls the Medical Text Indexer is invoked from this file

import json
import subprocess
import datetime
import os
import time
import config
import citations

# Creates a batch ID for .txt files to be used as input for MetaMap
def create_batch_id():
    # removes last 4 digits so ms is 2 s.f.
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-4] 

def format_results(results):
    for r in results:
        print(r)

def write_file(filename, results):
    batch = open(os.path.join('app/static', filename), 'a+') # default: unbuffered
    for i in range(0, len(results)):
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

class MetaMap():

    def run(self, results):
        batch_id = create_batch_id()
        filename = batch_id + '.txt'
        write_file(filename, results)

        mp = config.metamap_path

        # shell args for running MetaMapCaller.jar
        popen_args = ['java',
                      '-cp',
                      mp + '/classes:' +
                      mp + '/lib/skrAPI.jar:' +
                      mp + '/lib/commons-logging-1.1.1.jar:' +
                      mp + '/lib/httpclient-cache-4.1.1.jar:' +
                      mp + '/lib/httpcore-nio-4.1.jar:' +
                      mp + '/lib/httpclient-4.1.1.jar:' +
                      mp + '/lib/httpcore-4.1.jar:' +
                      mp + '/lib/httpmime-4.1.1.jar:' +
                      '.:' +
                      config.json_simple_path,
                      'MetaMapCaller',
                      '../../static/' + filename, 
                      config.un, # username for MetaMap
                      config.pwd, # password for MetaMap
                      config.email, # email for MetaMap
                      '-y'] 

        p = subprocess.Popen(popen_args, cwd='app/java/bin/', stdout=subprocess.PIPE)

        terms_list = list()

        while 1:
            term = p.stdout.readline()
            if not term and p.returncode is not None:
                break
            terms_list.append(term.decode('UTF-8'))
            p.poll()
        print("done %d" % p.returncode)

        format_results(terms_list)

# # testing..
if __name__ == "__main__":
    with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 
            'r') as datafile:
        mm = MetaMap()
        mm.run(json.load(datafile))
