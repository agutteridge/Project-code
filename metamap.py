# The Java .jar that calls the Medical Text Indexer is invoked from this file

import json, subprocess, datetime, os
import config, citations

# Creates a batch ID for .txt files to be used as input for MetaMap
def create_batch_id():
    # removes last 4 digits so ms is 2 s.f.
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-4] 

def run(results):
    batch_id = create_batch_id()
    filename = batch_id + '.txt'
    batch = open(filename, 'a+') # default: unbuffered
    for i in range(0, len(results)):
        ASCII_title = results[i]['MedlineCitation']['Article']['ArticleTitle'].encode('ascii', 
            errors='ignore').decode('UTF-8')
        ASCII_abstract = results[i]['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode('ascii', 
            errors='ignore').decode('UTF-8')
        batch.write('UI  - ' + batch_id + str(i) + '\n' + 
            'TI  - ' + ASCII_title + '\n' +
            'AB  - ' + ASCII_abstract + 
            '\n\n')

    print('done!')
    done = True

    mp = config.metamap_path

    # shell args for running MetaMapCaller.jar
    popen_args = ['java',
                  '-cp',
                  '../classes:' +
                  '../lib/skrAPI.jar:' +
                  '../lib/commons-logging-1.1.1.jar:' +
                  '../lib/httpclient-cache-4.1.1.jar:' +
                  '../lib/httpcore-nio-4.1.jar:' +
                  '../lib/httpclient-4.1.1.jar:' +
                  '../lib/httpcore-4.1.jar:' +
                  '../lib/httpmime-4.1.1.jar:' +
                  config.json_simple_path,
                  '-jar',
                  'MetaMapCaller.jar', 
                  '../../' + filename, 
                  config.email]

    for a in popen_args:
      print(a)      

    if done:
        p = subprocess.Popen(popen_args, cwd='./java/lib/', stdout=subprocess.PIPE)
        # JSON obj!
        output = p.stdout
        for line in output:
            print(line)

# testing..
if __name__ == "__main__":
    with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 
            'r') as datafile:
        run(json.load(datafile))