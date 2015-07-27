# The Java .jar that calls the Medical Text Indexer is invoked from this file

import json, subprocess, config

def main(results):
	filename = create_batch_id() + '.txt'
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

	if done:
		p = subprocess.Popen(['java', '-jar', 'MetaMapCaller.jar', filename, config.email], 
			cwd='./java/lib/', stdout=subprocess.PIPE)
		# JSON obj!
		print(p.stdout.readline().decode('UTF-8'))
