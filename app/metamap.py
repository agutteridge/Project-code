# The Java .jar that calls the Medical Text Indexer is invoked from this file

import subprocess
import datetime
import os
import time
import re

from app import config

# Creates a batch ID for .txt files to be used as input for MetaMap
def create_batch_id():
    # removes last 4 digits so ms is 2 s.f.
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-4] 

# Formats MetaMap results into a list of dicts, each with a PubMed ID and a list of 
# lists of concepts for that paper
def format_results(results, is_local):
    all_concepts = []
    for r in results:
        concept = r.split('|')
    
        if len(concept) > 1:
            pmid = concept[0]
            pmid_non_digits = re.sub('[\d]', '', pmid)
    
            # Local MetaMap and MetaMap via the Web API format results slightly differently
            if is_local:
                rest = [concept[3], concept[4]]
            else:
                rest = [concept[1], concept[2]]

            # filtering out empty strings and those with non-digits
            if len(pmid_non_digits) is 0 and len(pmid) is not 0:
                if not all_concepts or all_concepts[0]['PMID'] != pmid:
                    # keeping structure for umls.run
                    all_concepts = [{'PMID': pmid,
                        'concepts': []}] + all_concepts # prepend

                all_concepts[0]['concepts'].append(rest) # adds concept to list of concepts for that paper
    
    return all_concepts

# creates a text file for MetaMap to use as a source
def write_file(filename, results, is_local):
    batch = open(os.path.join('app/static', filename), 'w') # default: unbuffered
    for i in range(0, len(results)):
        # exclude ASCII characters to avoid MetaMap errors
        ASCII_title = results[i]['MedlineCitation']['Article']['ArticleTitle'].encode('ascii', 
            errors='ignore').decode('UTF-8')
        ASCII_abstract = ''

        # Local MetaMap returns too many concepts if abstract is included
        if not is_local:
            # not all papers have abstracts
            if 'Abstract' in results[i]['MedlineCitation']['Article']:
                ASCII_abstract = results[i]['MedlineCitation']['Article']['Abstract']['AbstractText'][0].encode('ascii', 
                    errors='ignore').decode('UTF-8')

        # Appending author keywords to abstract
        if 'KeywordList' in results[i]['MedlineCitation']:
            for kw in results[i]['MedlineCitation']['KeywordList']:
                for k in kw:
                    ASCII_abstract = ASCII_abstract + ' ' + str(k).encode('ascii', 
                        errors='ignore').decode('UTF-8') + '. '

        # Appending MeSH keywords to abstract
        if 'MeshHeadingList' in results[i]['MedlineCitation']:
            for mh in results[i]['MedlineCitation']['MeshHeadingList']:
                ASCII_abstract = ASCII_abstract + ' ' + mh['DescriptorName'].encode('ascii', 
                    errors='ignore').decode('UTF-8') + '. '

        batch.write('UI  - ' + results[i]['MedlineCitation']['PMID'] + '\n' + 
            'TI  - ' + ASCII_title + '\n' +
            'AB  - ' + ASCII_abstract + 
            '\n\n')

    batch.close()
    return True

def read_local_results():
    with open(os.path.join('app/static', 'mm_output.txt'), 'r') as datafile:
        txt = datafile.read()
        lines = txt.split('\n')
        
        datafile.close()
        return format_results(lines, True)

def run_local(results):
    print('in metamap.run (local)')

    if len(results) > 0:
        batch_id = create_batch_id()
        filename = batch_id + '.txt'
        write_file(filename, results, True)

        mp = config.mm_local_path

        # shell args for running MetaMapCaller.class
        popen_args = [mp + '/bin/metamap14',
                      '--strict_model',
                      '--prefer_multiple_concepts',
                      '--fielded_mmi_output',
                      '--prune',
                      '1',
                      config.proj_path + '/app/static/' + filename,
                      config.proj_path + '/app/static/mm_output.txt']
        
        p = subprocess.Popen(popen_args, stdout=subprocess.PIPE)

        while 1:
            term = p.stdout.readline()
            print(term)
            if not term and p.returncode is not None:
                break
            p.poll()
        print("done %d" % p.returncode)

        # Delete input file from app/static
        os.remove(os.path.join('app/static', filename))

        print('returning from metamap.run')
        return read_local_results()
    else:
        return []

def run(results):
    print('in metamap.run')

    if len(results) > 0:
        batch_id = create_batch_id()
        filename = batch_id + '.txt'
        write_file(filename, results, False)

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
            print(term)
            if not term and p.returncode is not None:
                break
            terms_list.append(term.decode('UTF-8'))
            p.poll()
        print("done %d" % p.returncode)

        # Delete input file from app/static
        os.remove(os.path.join('./app/static', filename))

        print('returning from metamap.run')
        end(now)
        return format_results(terms_list, False)
    else:
        return []
