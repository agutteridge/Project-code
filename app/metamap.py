# The Java .jar that calls the Medical Text Indexer is invoked from this file

import json, subprocess, datetime, os, time
import config, citations

# Creates a batch ID for .txt files to be used as input for MetaMap
def create_batch_id():
    # removes last 4 digits so ms is 2 s.f.
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-4] 

def run(results):
    batch_id = create_batch_id()
    filename = batch_id + '.txt'
    batch = open(os.path.join('./java/src', filename), 'a+') # default: unbuffered
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
                  filename, 
                  config.un, # username for MetaMap
                  config.pwd,
                  config.email] # password for MetaMap

    if done:
        p = subprocess.Popen(popen_args, cwd='./java/bin/', stdout=subprocess.PIPE)
        # time.sleep(1)
        # print(config.email)
        # # time.sleep(10)
        # print('oh hey password')
        # print('ohheypasswrod')
        # p.stdin.write(bytes(config.email + '\n', 'UTF-8'))
        # p.communicate(input=bytes(config.pwd + '\n', 'UTF-8'))
        # p.communicate(input=bytes(config.email + '\n', 'UTF-8'))
        # p.stdin.write(bytes(config.email + '\n', 'UTF-8'))

        # p.stdin.write('\t' + config.email + '\n')
        # p.communicate('\t' + config.pwd + '\n')
        terms_list = list()

        terms_output = p.stdout.readline()
        
        while terms_output is not 'END\n':
            if terms_output is not '':
                print(terms_output)
                terms_list.extend(terms_output)
                terms_output = p.stdout.readline()

        # for line in enumerate(lines_iterator):
        #     print('line number is: ' + str(line[0]))
        #     if line[0] is 0:
        #         p.stdin.write(config.email)
        #         p.communicate
        #     elif line[0] is 1:
        #         p.stdin.write(config.pwd)
        #     else:
        #         # MetaMap data
        #         print(line[1])

# testing..
if __name__ == "__main__":
    with open(os.path.join('../tests/resources', 'eFetch_sample.json'), 
            'r') as datafile:
        run(json.load(datafile))