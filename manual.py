import os
import json

from app import citations, geocode

def fake_json(filename):
    with open(os.path.join('./tests/resources', filename), 'r') as datafile:
        obj = json.load(datafile)
        datafile.close()
        return obj

def format_address_0(address):
    return address

def geocode_run(results):
    with open(os.path.join('./tests/resources', 'format_address_0_results.txt'), 'a') as datafile:
        input_num = 0
        output_num = 0

        for paper in results:
            pmid = str(paper['MedlineCitation']['PMID'])
            datafile.write('PMID: ' + pmid + '\n')

            if 'AuthorList' in paper['MedlineCitation']['Article']: 
                author_list = paper['MedlineCitation']['Article']['AuthorList']

                for author in author_list:
                    if 'LastName' in author:
                        datafile.write('\tAuthor: ' + author['LastName'] + '\n')

                    if 'CollectiveName' in author:
                        datafile.write('\tCollective: ' + author['CollectiveName'] + '\n')

                    for place in author['AffiliationInfo']:
                        individual_addresses = place['Affiliation'].split(';')

                        for f in individual_addresses:
                            datafile.write('\t\tInput address: ' + f + '\n')
                            input_num += 1
                            formatted_address = format_address_0(f) # change format_address
                            place = geocode.get_location(formatted_address)
        
                            if place:
                                datafile.write('\t\tOutput address: ' + place['name'] + '\n')
                                output_num += 1
                            else:
                                datafile.write('\t\tOutput address: NONE' + '\n')
            else:
                datafile.write('\tNO AUTHORS LISTED.')

        datafile.write('\nSUCCESS RATE: ' + str(output_num / input_num * 100) + '%')
    datafile.close()


# Writes first 100 papers from PubMed retrieved with query 'biology' to JSON file 
def init():
    with open(os.path.join('./tests/resources', 'eFetch_random.json'), 'a') as datafile:
        print('opened')
        count = 0
        json_list = []

        while count < 100:
            print('loop')
            id_list = citations._search('biology')['IdList']
            fetch_results = citations.fetch_details(id_list)

            for fr in fetch_results:
                count += len(fetch_results)
                json_list.append(fr)

        datafile.write(json.dumps(json_list))
    datafile.close()

if __name__ == "__main__":
    geocode_run(fake_json('eFetch_random.json'))
    # print('no method chosen')
